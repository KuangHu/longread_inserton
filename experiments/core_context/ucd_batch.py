"""Batch run of the U-insertion-D comparison over an existing candidate catalogue, with two checks per pair.

Scope: a full comparison of the pairs of the existing candidate catalogue (here every match_both_lt2500 pair of the
boundary review, b1-2-3). It is not a complete discovery over every FASTA. Pairs are taken as listed; nothing is
filtered by old read-check, gap-caller or S6 results (s6_grouping is carried as a label only).

Per pair and interval kind (cds_core, candidate_insertion; always output separately):
  screen         ucd_compare.compare (v0.1, frozen), on the sequences the screen extracts itself. Kept as is.
  source_check   Check 1 (ucd_checks.check1_side). Both sides must pass:
                 - the source FASTA sha256 equals the verified manifest;
                 - the contig exists and the interval is in bounds;
                 - U/X/D read again by byte offset equal the screen's, base by base, with the minus strand
                   reverse-complemented.
  ucd_check      Check 2 (ucd_checks.check_ucd) on the Check-1 sequences: collinear insertion chain, frozen context rule.
  final_status   pass = source_check pass and ucd_check ucd_supported; otherwise not_passed, with failure_reason.

Commands
  freeze   build the input pair table from the catalogue + boundary / occurrence shards; verify each source FASTA against
           the manifest; split into chunks; record sha256 of code, inputs and parameters (freeze.json).
  chunk    process one chunk (SLURM array index). A finished chunk (done.json matching the freeze) is skipped, so a
           resubmitted array resumes where it stopped. Code changed since freeze -> refuse to run.
  merge    per-ID completeness (every pair_id x interval kind exactly once), then final tables, counts, examples.
  cleanup  after a complete merge: remove the chunk directories (rebuildable). Keeps the source table, the final
           sequences and the alignment evidence.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import socket
import sys
import tempfile
import time
from collections import Counter, OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import ucd_checks as ck
import ucd_compare as uc

VERSION = 'ucd_batch_v0.1'
HERE = Path(__file__).resolve().parent
CODE = ('ucd_compare.py', 'ucd_checks.py', 'ucd_batch.py')
KINDS = uc.KINDS
SIDE = ('occurrence', 'sample', 'assembly', 'contig', 'start0', 'end0', 'strand', 'ext_up', 'ext_down', 'elem_lo',
        'elem_hi', 'boundary_up', 'boundary_down', 'estimated_span', 'input_error')
LABELS = ('pair_id', 'structure', 'core_key', 'pair_scope', 'copy_class', 's6_grouping', 'core_match_class',
          'core_dna_identity', 'core_dna_coverage_min', 'protein_identity_min')
INPUT_COLS = list(LABELS) + ['chunk'] + ['%s_%s' % (s, f) for s in 'ab' for f in SIDE]
SCREEN = ('verdict', 'insertion_class', 'insertion_identity', 'insertion_cov_a', 'insertion_cov_b', 'b_flipped',
          'UU_class', 'UU_cov_a', 'UU_cov_b', 'DD_class', 'DD_cov_a', 'DD_cov_b', 'reasons')
PROV = ('occurrence', 'sample', 'assembly', 'fasta_sha256', 'contig', 'contig_len', 'strand', 'x_interval',
        'u_interval', 'd_interval', 'edge', 'x_sha', 'u_sha', 'd_sha', 'boundary_up', 'boundary_down',
        'estimated_span')
CHK = (['chk_b_flipped', 'chk_insertion_class', 'chk_insertion_identity', 'chk_insertion_cov_a', 'chk_insertion_cov_b',
        'chk_insertion_chain_hsps', 'chk_insertion_orientation', 'chk_insertion_union_cov_a',
        'chk_insertion_union_cov_b', 'chk_U_A_usable', 'chk_D_A_usable', 'chk_U_B_usable', 'chk_D_B_usable'] +
       ['chk_%s_%s' % (c, f) for c, *_ in ck.PAIRS8 for f in ('class', 'identity', 'cov_a', 'cov_b')] +
       ['chk_upstream', 'chk_downstream', 'chk_reasons'])
ROW_COLS = (['pair_key', 'pair_id', 'interval_kind'] + list(LABELS[1:]) +
            ['%s_%s' % (s, f) for s in 'ab' for f in PROV] + ['x_len_a', 'x_len_b'] +
            ['screen_' + k for k in SCREEN] + ['source_check', 'source_reasons', 'ucd_check'] + CHK +
            ['final_status', 'failure_reason', 'screen_check_agree'])
HSP_COLS = ['pair_key', 'q', 's', 'strand', 'pident', 'length', 'nident', 'qs', 'qe', 'ss', 'se', 'qlen', 'slen',
            'evalue', 'bitscore', 'in_insertion_chain']


def sha_file(p):
    return uc.sha256(p)


def code_hashes():
    return {c: sha_file(HERE / c) for c in CODE}


def pair_id(structure, core, scope, oa, ob):
    return 'M' + hashlib.sha1('|'.join((structure, core, scope, oa, ob)).encode()).hexdigest()[:12]


# ---------------------------------------------------------------- freeze

def cmd_freeze(a):
    out, work = Path(a.out), Path(a.work)
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    cat = uc.read_tsv(a.catalogue)
    keys = {}
    for r in cat:
        oa, ob = sorted((r['occurrence_a'], r['occurrence_b']))
        k = (r['structure'], r['core_key'], r['pair_scope'], oa, ob)
        if k in keys:
            sys.exit('duplicate catalogue pair %s' % (k,))
        keys[k] = r
    found, need = {}, set()
    for i in range(64):
        for r in uc.read_tsv(Path(a.review) / 'boundary' / a.tag / ('shard_%03d' % i) / 'pairs.tsv.gz'):
            oa, ob = r['occurrence_a'], r['occurrence_b']
            sw = oa > ob
            k = (r['structure'], r['core_key'], r['pair_scope']) + ((ob, oa) if sw else (oa, ob))
            if k in keys:
                if k in found:
                    sys.exit('pair found twice in boundary shards %s' % (k,))
                found[k] = (r, sw)
                need |= {oa, ob}
    occ = {}
    for i in range(64):
        for r in uc.read_tsv(Path(a.review) / 's6' / a.tag / ('shard_%03d' % i) / 'occurrences.tsv.gz'):
            if r['occurrence_id'] in need:
                occ[(r['structure'], r['occurrence_id'])] = r
    bnd = {}
    for scope in ('same_assembly', 'cross_assembly'):
        for r in uc.read_tsv(uc.REPO / ('boundary_occurrences_%s_%s.tsv.gz' % (scope, a.tag))):
            if r['occurrence_id'] in need:
                bnd[(scope, r['structure'], r['occurrence_id'])] = r
    man = uc.manifest(a.manifest)
    rows = []
    for k, c in keys.items():
        structure, core, scope, oa, ob = k
        x = dict(pair_id=pair_id(*k), structure=structure, core_key=core, pair_scope=scope,
                 **{f: c.get(f, '') for f in LABELS[4:]})
        f, sw = found.get(k, (None, False))
        for s, o in (('a', oa), ('b', ob)):
            e = []
            if f is None:
                e.append('pair_not_in_boundary_review')
            ro = occ.get((structure, o))
            if ro is None:
                e.append('occurrence_coordinates_missing')
            d = dict(occurrence=o, sample=o.split('|')[1], input_error=','.join(e))
            if ro:
                d.update(sample=ro['sample'], contig=ro['contig'], start0=ro['start0'], end0=ro['end0'],
                         strand=ro['strand'])
            d['assembly'] = man.get(d['sample'], {}).get('assembly', '')
            if not d['assembly']:
                d['input_error'] = ','.join(filter(None, (d['input_error'], 'sample_not_in_manifest')))
            if f is not None and ro:
                src = ('b' if s == 'a' else 'a') if sw else s
                up, down = f[src + '_up'], f[src + '_down']
                d.update(ext_up=up, ext_down=down)
                if up != '' and down != '':
                    st, en = int(ro['start0']), int(ro['end0'])
                    lo, hi = (st - int(up), en + int(down)) if ro['strand'] == '+' else (st - int(down), en + int(up))
                    d.update(elem_lo=lo, elem_hi=hi)
            b = bnd.get((scope, structure, o), {})
            d.update(boundary_up=b.get('up_status', ''), boundary_down=b.get('down_status', ''),
                     estimated_span=b.get('estimated_span', ''))
            x.update({'%s_%s' % (s, q): v for q, v in d.items()})
        rows.append(x)
    rows.sort(key=lambda r: (r['a_sample'], r['b_sample'], r['pair_id']))
    for i, r in enumerate(rows):
        r['chunk'] = i // a.chunk_size
    n_chunks = (len(rows) + a.chunk_size - 1) // a.chunk_size
    uc.write_tsv(out / 'input_pairs.tsv.gz', rows, INPUT_COLS)
    # sources: every assembly used, verified against the manifest sha256 now (compute nodes re-verify on load)
    samples = sorted({r[s + '_sample'] for r in rows for s in 'ab'})
    src = []
    for smp in samples:
        m = man.get(smp)
        if m is None:
            src.append(dict(sample=smp, status='not_in_manifest'))
            continue
        p = Path(m['fasta'])
        h = sha_file(p) if p.exists() else ''
        src.append(dict(sample=smp, assembly=m['assembly'], fasta=m['fasta'], manifest_sha256=m['fna_sha256'],
                        observed_sha256=h, bytes=p.stat().st_size if p.exists() else '',
                        status='verified' if h == m['fna_sha256'] else ('missing' if not h else 'sha256_mismatch')))
    uc.write_tsv(out / 'sources.tsv.gz', src, ['sample', 'assembly', 'fasta', 'manifest_sha256', 'observed_sha256',
                                               'bytes', 'status'])
    cfg = dict(uc.CONFIG)
    frz = dict(version=VERSION, screen_version=uc.VERSION, check_version=ck.CHECK_VERSION, config=cfg,
               check_config=ck.CHECK_CONFIG, coordinates=ck.COORDS, code_sha256=code_hashes(),
               scope='full comparison of the existing candidate catalogue (match_both_lt2500, %s); not a complete '
                     'discovery over all FASTA' % a.tag,
               inputs={str(p): sha_file(p) for p in (Path(a.catalogue), Path(a.manifest))},
               input_pairs=len(rows), input_pairs_sha256=sha_file(out / 'input_pairs.tsv.gz'),
               sources_sha256=sha_file(out / 'sources.tsv.gz'), chunk_size=a.chunk_size, chunks=n_chunks,
               input_errors=dict(Counter(e for r in rows for s in 'ab' for e in r[s + '_input_error'].split(',') if e)),
               sources=dict(Counter(s['status'] for s in src)), work=str(work),
               frozen_at=time.strftime('%Y-%m-%d %H:%M:%S'))
    (out / 'freeze.json').write_text(json.dumps(frz, indent=1) + '\n')
    print(json.dumps({k: v for k, v in frz.items() if k != 'config'}, indent=1))


# ---------------------------------------------------------------- chunk

_W = {}


def _init(man, cfg, ccfg):
    _W.update(man=man, cfg=cfg, ccfg=ccfg, genomes=OrderedDict(), sources=OrderedDict())


def _lru(cache, key, make, size=48):
    if key in cache:
        cache.move_to_end(key)
        return cache[key]
    v = cache[key] = make()
    if len(cache) > size:
        cache.popitem(last=False)
    return v


def screen_genome(sample):
    return _lru(_W['genomes'], sample, lambda: uc.read_fasta(_W['man'][sample]['fasta']))


def source(sample):
    m = _W['man'][sample]
    return _lru(_W['sources'], sample, lambda: ck.Source(m['fasta'], m['fna_sha256']))


def interval(r, s, kind):
    if r[s + '_input_error']:
        return 'input_error:' + r[s + '_input_error']
    if kind == 'cds_core':
        return int(r[s + '_start0']), int(r[s + '_end0'])
    if r[s + '_elem_lo'] == '' or r[s + '_boundary_up'] != 'resolved' or r[s + '_boundary_down'] != 'resolved':
        return None
    return int(r[s + '_elem_lo']), int(r[s + '_elem_hi'])


def work_pair(r):
    cfg, ccfg, L = _W['cfg'], _W['ccfg'], _W['cfg']['context']
    rows, seqs, hsps = [], {}, []
    t0 = sum(os.times()[:4])  # own + child (blastn) CPU
    tw = time.time()
    tmp = tempfile.mkdtemp(prefix='ucdb_', dir=os.environ.get('TMPDIR'))
    try:
        for kind in KINDS:
            key = '%s|%s' % (r['pair_id'], kind)
            row = dict({k: r[k] for k in LABELS}, pair_key=key, interval_kind=kind)
            for s in 'ab':
                row.update({'%s_%s' % (s, f): r[s + '_' + f] for f in ('occurrence', 'sample', 'assembly', 'contig',
                                                                         'strand', 'boundary_up', 'boundary_down',
                                                                         'estimated_span')})
            ivs = {s: interval(r, s, kind) for s in 'ab'}
            bad = [s for s in 'ab' if not isinstance(ivs[s], tuple)]
            if bad:
                why = ';'.join('%s:%s' % (s, ivs[s] or 'interval_not_available') for s in bad)
                row.update(source_check='not_run', source_reasons=why, ucd_check='not_run',
                           final_status='not_passed', failure_reason=why)
                rows.append(row)
                continue
            # screen (v0.1, its own extraction)
            S = {}
            for s in 'ab':
                g = screen_genome(r[s + '_sample'])
                cs = g.get(r[s + '_contig'])
                if cs is None:
                    S[s] = None
                    continue
                S[s] = dict(uc.sides_of(cs, ivs[s][0], ivs[s][1], r[s + '_strand'], L), id=r[s + '_occurrence'],
                            sample=r[s + '_sample'], contig=r[s + '_contig'])
            if S['a'] and S['b']:
                sc = uc.compare(S['a'], S['b'], cfg, tmp)
                row.update({'screen_' + k: sc[k] for k in SCREEN})
            else:
                row['screen_verdict'] = 'screen_error:contig_missing'
            # Check 1
            C, why = {}, []
            for s in 'ab':
                c = ck.check1_side(source(r[s + '_sample']), r[s + '_contig'], ivs[s][0], ivs[s][1],
                                   r[s + '_strand'], L, S[s])
                C[s] = c
                row.update({'%s_fasta_sha256' % s: getattr(source(r[s + '_sample']), 'sha256', ''),
                            '%s_contig_len' % s: c['contig_len'], '%s_edge' % s: c['edge']})
                if 'ivs' in c:
                    row.update({'%s_%s_interval' % (s, p.lower()): '%d-%d' % c['ivs'][p] for p in 'XUD'})
                    row.update({'%s_%s_sha' % (s, p.lower()): c['sha'][p] for p in 'XUD'})
                why += ['%s:%s' % (s, x) for x in c['reasons']]
            row['source_check'] = 'pass' if not why else 'fail'
            row['source_reasons'] = ';'.join(why)
            if why:
                row.update(ucd_check='not_run:source_check_failed', final_status='not_passed',
                           failure_reason='source_check:' + ';'.join(why))
                rows.append(row)
                continue
            # Check 2 on the verified sequences
            O = {s: dict(C[s]['seqs'], id=r[s + '_occurrence'], sample=r[s + '_sample'], contig=r[s + '_contig'],
                         ivs=C[s]['ivs']) for s in 'ab'}
            res, h = ck.check_ucd(O['a'], O['b'], cfg, ccfg, tmp)
            assert res['chk_occurrence_a'] == r['a_occurrence']  # input is in canonical order
            row.update({k: res[k] for k in CHK + ['ucd_check']})
            row['x_len_a'], row['x_len_b'] = len(O['a']['X']), len(O['b']['X'])
            ok = res['ucd_check'] == 'ucd_supported'
            row['final_status'] = 'pass' if ok else 'not_passed'
            row['failure_reason'] = '' if ok else 'ucd_check:%s (%s)' % (res['ucd_check'], res['chk_reasons'])
            row['screen_check_agree'] = row.get('screen_verdict') == res['ucd_check']
            rows.append(row)
            hsps += [dict(x, pair_key=key) for x in h]
            for s in 'ab':
                c = C[s]
                sid = '%s|%s|%d-%d' % (r[s + '_occurrence'], kind, *c['ivs']['X'])
                for p in 'UXD':
                    lo, hi = c['ivs'][p]
                    seqs['%s|%s' % (sid, p)] = ('%s|%s assembly=%s contig=%s interval=%d-%d strand=%s part=%s '
                                                'fasta_sha256=%s seq_sha256_16=%s' % (
                                                    sid, p, r[s + '_assembly'], r[s + '_contig'], lo, hi,
                                                    r[s + '_strand'], p, row[s + '_fasta_sha256'], c['sha'][p]),
                                                c['seqs'][p])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return rows, seqs, hsps, sum(os.times()[:4]) - t0, time.time() - tw


def cmd_chunk(a):
    out, work = Path(a.out), Path(a.work)
    frz = json.loads((out / 'freeze.json').read_text())
    if code_hashes() != frz['code_sha256']:
        sys.exit('code differs from the frozen code (freeze.json); refusing to run')
    idx = a.index if a.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID'])
    d = work / 'chunks' / ('c%04d' % idx)
    done = d / 'done.json'
    if done.exists():
        j = json.loads(done.read_text())
        if j.get('input_pairs_sha256') == frz['input_pairs_sha256'] and \
                sha_file(d / 'rows.tsv.gz') == j.get('rows_sha256'):
            print('chunk %d already complete; skipped' % idx)
            return
    pairs = [r for r in uc.read_tsv(out / 'input_pairs.tsv.gz') if int(r['chunk']) == idx]
    man = uc.manifest(a.manifest)
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows, seqs, hsps, cpu, wall = [], {}, [], [], []
    with ProcessPoolExecutor(a.threads, initializer=_init, initargs=(man, frz['config'], frz['check_config'])) as ex:
        for rs, sq, hs, c, w in ex.map(work_pair, pairs, chunksize=4):
            rows += rs
            seqs.update(sq)
            hsps += hs
            cpu.append(c)
            wall.append(w)
    tmpd = d / 'partial'
    tmpd.mkdir(exist_ok=True)
    uc.write_tsv(tmpd / 'rows.tsv.gz', rows, ROW_COLS)
    with gzip.open(tmpd / 'sequences.fa.gz', 'wt') as f:
        for k in sorted(seqs):
            f.write('>%s\n%s\n' % seqs[k])
    uc.write_tsv(tmpd / 'hsps.tsv.gz', hsps, HSP_COLS)
    for n in ('rows.tsv.gz', 'sequences.fa.gz', 'hsps.tsv.gz'):
        os.replace(tmpd / n, d / n)
    tmpd.rmdir()
    j = dict(chunk=idx, pairs=len(pairs), rows=len(rows), pair_ids_sha256=hashlib.sha256(
        '\n'.join(sorted(r['pair_id'] for r in pairs)).encode()).hexdigest(),
        rows_sha256=sha_file(d / 'rows.tsv.gz'), input_pairs_sha256=frz['input_pairs_sha256'],
        wall_s=round(time.time() - t0, 1), threads=a.threads, worker_wall_s=round(sum(wall), 1),
        child_process_cpu_s=round(sum(cpu), 1), host=socket.gethostname(),
        slurm=dict(job=os.environ.get('SLURM_ARRAY_JOB_ID', os.environ.get('SLURM_JOB_ID', '')),
                   task=os.environ.get('SLURM_ARRAY_TASK_ID', ''), restart=os.environ.get('SLURM_RESTART_COUNT', '')),
        finished=time.strftime('%Y-%m-%d %H:%M:%S'))
    done.write_text(json.dumps(j, indent=1) + '\n')
    print(json.dumps(j))


# ---------------------------------------------------------------- merge

def counts(rows):
    occ = {o for r in rows for o in (r['a_occurrence'], r['b_occurrence'])}
    return dict(pairs=len(rows), occurrences=len(occ),
                samples=len({r['a_sample'] for r in rows} | {r['b_sample'] for r in rows}),
                cores=len({(r['structure'], r['core_key']) for r in rows}))


def cmd_merge(a):
    out, work = Path(a.out), Path(a.work)
    frz = json.loads((out / 'freeze.json').read_text())
    inp = uc.read_tsv(out / 'input_pairs.tsv.gz')
    expect = {(r['pair_id'], k) for r in inp for k in KINDS}
    rows, seen, missing_chunks, done = [], Counter(), [], []
    for i in range(frz['chunks']):
        d = work / 'chunks' / ('c%04d' % i)
        if not (d / 'done.json').exists():
            missing_chunks.append(i)
            continue
        j = json.loads((d / 'done.json').read_text())
        if sha_file(d / 'rows.tsv.gz') != j['rows_sha256']:
            missing_chunks.append(i)
            continue
        done.append(j)
        for r in uc.read_tsv(d / 'rows.tsv.gz'):
            seen[(r['pair_id'], r['interval_kind'])] += 1
            rows.append(r)
    missing = sorted(expect - set(seen))
    dup = sorted(k for k, v in seen.items() if v > 1)
    extra = sorted(set(seen) - expect)
    comp = dict(expected=len(expect), observed=len(seen), missing=len(missing), duplicated=len(dup), unexpected=len(extra),
                missing_chunks=missing_chunks, missing_examples=missing[:20], duplicated_examples=dup[:20],
                complete=not (missing or dup or extra or missing_chunks))
    (out / 'completeness.json').write_text(json.dumps(comp, indent=1) + '\n')
    if not comp['complete']:
        print(json.dumps(comp, indent=1))
        sys.exit('incomplete: rerun the missing chunks (the array resumes) before merging')
    res = {}
    for kind in KINDS:
        rk = sorted((r for r in rows if r['interval_kind'] == kind), key=lambda r: r['pair_id'])
        uc.write_tsv(work / ('ucd_all_rows_%s.tsv.gz' % kind), rk, ROW_COLS)
        uc.write_tsv(out / ('ucd_final_%s.tsv.gz' % kind), [r for r in rk if r['final_status'] == 'pass'], ROW_COLS)
        by = defaultdict(list)
        for r in rk:
            by[r['ucd_check']].append(r)
        res[kind] = dict(
            rows=len(rk),
            source_check=dict(Counter(r['source_check'] for r in rk)),
            source_failure_reasons=dict(Counter(x.split(':', 1)[1].split('@')[0] for r in rk if r['source_check'] == 'fail'
                                                for x in r['source_reasons'].split(';'))),
            not_run_reasons=dict(Counter(x.split(':', 1)[1] for r in rk if r['source_check'] == 'not_run'
                                         for x in r['source_reasons'].split(';'))),
            ucd_check={v: counts(x) for v, x in sorted(by.items())},
            final_pass=counts([r for r in rk if r['final_status'] == 'pass']),
            final_pass_by_scope_structure={'%s|%s' % k: counts(x) for k, x in sorted(_group(
                [r for r in rk if r['final_status'] == 'pass'], ('pair_scope', 'structure')).items())},
            screen_vs_check=dict(sorted(Counter('%s -> %s' % (r.get('screen_verdict') or '-', r['ucd_check'])
                                                for r in rk if r['source_check'] == 'pass').items())),
            screen_check_disagree=sum(r['screen_check_agree'] == 'False' for r in rk))
    # evidence: final sequences (deduplicated by occurrence interval) and Check-2 HSPs
    seq_p, hsp_p = work / 'evidence_sequences.fa.gz', work / 'evidence_hsps.tsv.gz'
    seen_s = set()
    with gzip.open(seq_p, 'wt') as fo:
        for i in range(frz['chunks']):
            with gzip.open(work / 'chunks' / ('c%04d' % i) / 'sequences.fa.gz', 'rt') as f:
                for h, s in zip(f, f):
                    k = h.split()[0]
                    if k not in seen_s:
                        seen_s.add(k)
                        fo.write(h + s)
    with gzip.open(hsp_p, 'wt') as fo:
        for i in range(frz['chunks']):
            with gzip.open(work / 'chunks' / ('c%04d' % i) / 'hsps.tsv.gz', 'rt') as f:
                hdr = f.readline()
                if i == 0:
                    fo.write(hdr)
                shutil.copyfileobj(f, fo)
    cpu = sum(j['child_process_cpu_s'] for j in done)
    summ = dict(
        version=VERSION, freeze=dict((k, frz[k]) for k in ('code_sha256', 'config', 'check_config', 'coordinates',
                                                           'scope', 'input_pairs_sha256', 'sources_sha256',
                                                           'input_errors', 'sources')),
        completeness=comp, by_kind=res,
        timing=dict(chunks=len(done), sum_chunk_wall_s=round(sum(j['wall_s'] for j in done), 1),
                    worker_cpu_s=round(cpu, 1), worker_wall_s=round(sum(j['worker_wall_s'] for j in done), 1),
                    worker_wall_s_per_pair=round(sum(j['worker_wall_s'] for j in done) / max(1, len(inp)), 3),
                    first_finished=min(j['finished'] for j in done), last_finished=max(j['finished'] for j in done),
                    hosts=dict(Counter(j['host'] for j in done))),
        files={p.name: dict(path=str(p), bytes=p.stat().st_size, sha256=sha_file(p)) for p in
               [work / ('ucd_all_rows_%s.tsv.gz' % k) for k in KINDS] + [seq_p, hsp_p] +
               [out / ('ucd_final_%s.tsv.gz' % k) for k in KINDS]})
    (out / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(dict(completeness=comp, by_kind={k: dict(final_pass=v['final_pass'], source_check=v['source_check'])
                                                      for k, v in res.items()}, timing=summ['timing']), indent=1))


def _group(rows, cols):
    g = defaultdict(list)
    for r in rows:
        g[tuple(r[c] for c in cols)].append(r)
    return g


# ---------------------------------------------------------------- cleanup

def cmd_cleanup(a):
    out, work = Path(a.out), Path(a.work)
    comp = json.loads((out / 'completeness.json').read_text())
    summ = json.loads((out / 'summary.json').read_text())
    if not comp.get('complete'):
        sys.exit('merge not complete; nothing removed')
    for name, f in summ['files'].items():
        if sha_file(f['path']) != f['sha256']:
            sys.exit('merged file %s does not match summary.json; nothing removed' % name)
    d = work / 'chunks'
    freed = sum(p.stat().st_size for p in d.rglob('*') if p.is_file())
    if a.apply:
        shutil.rmtree(d)
    rec = dict(removed=str(d), bytes=freed, applied=a.apply, kept=sorted(summ['files']) + ['sources.tsv.gz'],
               at=time.strftime('%Y-%m-%d %H:%M:%S'))
    (out / 'cleanup.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps(rec, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest='cmd', required=True)
    for name in ('freeze', 'chunk', 'merge', 'cleanup'):
        x = sp.add_parser(name)
        x.add_argument('--out', required=True)
        x.add_argument('--work', required=True)
        x.add_argument('--manifest', default=str(uc.REPO / 'manifest.tsv'))
        if name == 'freeze':
            x.add_argument('--catalogue', default=str(uc.REPO / 'match_both_lt2500_pairs_b1-2-3.tsv.gz'))
            x.add_argument('--review', default=str(uc.SCRATCH / 'scale_comamonadaceae_v0.1' / 'review'))
            x.add_argument('--tag', default='b1-2-3')
            x.add_argument('--chunk-size', type=int, default=2000)
        if name == 'chunk':
            x.add_argument('--index', type=int)
            x.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
        if name == 'cleanup':
            x.add_argument('--apply', action='store_true')
    a = p.parse_args()
    dict(freeze=cmd_freeze, chunk=cmd_chunk, merge=cmd_merge, cleanup=cmd_cleanup)[a.cmd](a)


if __name__ == '__main__':
    main()
