"""Batch run of the U-insertion-D comparison over an existing candidate catalogue, with two checks per pair.

Scope: a full comparison of the pairs of an existing candidate catalogue. It is not a complete discovery over every
FASTA. Pairs are taken as listed; nothing is filtered by old read-check, gap-caller or S6 results (s6_grouping is
carried as a label only).

Per pair and interval kind (cds_core, candidate_insertion; always output separately):
  screen         ucd_compare.compare (frozen rules), on the sequences the screen extracts itself. Kept as is.
  source_check   Check 1 (ucd_checks.check1_side). Both sides must pass:
                 - the source FASTA sha256 equals the verified manifest;
                 - the contig exists and the interval is in bounds;
                 - U/X/D read again by byte offset equal the screen's, base by base, with the minus strand
                   reverse-complemented.
  ucd_check      Check 2 (ucd_checks.check_ucd) on the Check-1 sequences: collinear insertion chain, frozen context rule.
  final_status   pass = source_check pass and ucd_check ucd_supported; otherwise not_passed, with failure_reason.
  A candidate_insertion interval exists only where the boundary review resolved both ends of both occurrences. Where
  it does not, the row says which end is missing; no boundary is inferred here.

Commands
  freeze   build the input table and split it into chunks; record sha256 of code, inputs and parameters (freeze.json).
           --catalogue match FILE     the pairs listed in a catalogue file (as v0.1)
           --catalogue boundary       every core pair of the boundary review shards; writes catalogue_ids.txt.gz
           --exclude-input FILE       drop the pair IDs of an earlier input table (ID difference, counted)
           --from-input FILE --ids F  reuse rows of an earlier frozen input table for the listed IDs (reruns)
  chunk    process one chunk (SLURM array index). A finished chunk (done.json matching the freeze) is skipped, so a
           resubmitted array resumes. Code changed since freeze -> refuse to run.
  merge    per-ID completeness of this run (every pair_id x interval kind exactly once); then rows per kind and
           evidence, streamed. Detailed tables and counts: ucd_combine.py.
  cleanup  after a complete merge: remove chunk directories and chunk inputs (rebuildable from freeze); writes the
           list of removed paths first.
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
from collections import Counter, OrderedDict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import ucd_checks as ck
import ucd_compare as uc

VERSION = 'ucd_batch_v0.4'  # v0.4: --input-table (pre-built pair table, e.g. the DNA entry) and --kinds
# The rule that turns input rows into intervals and verdicts into final_status (interval(), work_pair()). Unchanged
# since v0.2; bump it if that logic changes, so earlier rows stop being reusable.
DECISION_VERSION = 'decision_v0.2'
HERE = Path(__file__).resolve().parent
CODE = ('ucd_compare.py', 'ucd_checks.py', 'ucd_batch.py')
KINDS = uc.KINDS  # interval kinds of this run; a run may restrict them (freeze --kinds, recorded in freeze.json)


def set_kinds(kinds):
    global KINDS
    kinds = tuple(kinds)
    assert kinds and all(k in uc.KINDS for k in kinds), kinds
    KINDS = kinds
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
            ['final_status', 'failure_reason', 'screen_check_agree', 'code_version'])
HSP_COLS = ['pair_key', 'q', 's', 'strand', 'pident', 'length', 'nident', 'qs', 'qe', 'ss', 'se', 'qlen', 'slen',
            'evalue', 'bitscore', 'in_insertion_chain']


def sha_file(p):
    return uc.sha256(p)


def code_hashes():
    return {c: sha_file(HERE / c) for c in CODE}


def pair_id(structure, core, scope, oa, ob):
    return 'M' + hashlib.sha1('|'.join((structure, core, scope, oa, ob)).encode()).hexdigest()[:12]


def fingerprint(r, kind, man, cfg, ccfg, code):
    """Input fingerprint of one pair and interval kind: source files, intervals / reasons and strands of both sides,
    parameters, and the screen / check code. Equal fingerprints mean the earlier row can be reused as it is."""
    parts = [kind, DECISION_VERSION, json.dumps(cfg, sort_keys=True), json.dumps(ccfg, sort_keys=True),
             code['ucd_compare.py'], code['ucd_checks.py']]
    for s in 'ab':
        smp = r[s + '_sample']
        parts += [smp, r[s + '_assembly'], man.get(smp, {}).get('fna_sha256', 'no_manifest_entry'), r[s + '_contig'],
                  r[s + '_strand'], str(interval(r, s, kind))]
    return hashlib.sha256('\x1f'.join(map(str, parts)).encode()).hexdigest()[:24]


def ids_of(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        first = f.readline().rstrip('\n').split('\t')
        if 'pair_id' in first:
            i = first.index('pair_id')
            return {line.split('\t')[i].rstrip('\n') for line in f}
        return {first[0]} | {line.strip() for line in f if line.strip()}


# ---------------------------------------------------------------- freeze

def build_from_catalogue(a, man):
    """Input rows (lists in INPUT_COLS order) for the catalogue pairs, from the boundary and occurrence shards."""
    keys = None
    if a.catalogue != 'boundary':
        keys = {}
        for r in uc.read_tsv(a.catalogue):
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
            if keys is not None and k not in keys:
                continue
            if k in found:
                sys.exit('pair found twice in boundary shards %s' % (k,))
            lab = keys[k] if keys is not None else r
            found[k] = ([lab.get(f, '') for f in LABELS[4:]], sw, r['a_up'], r['a_down'], r['b_up'], r['b_down'])
            need |= {oa, ob}
    occ = {}
    for i in range(64):
        for r in uc.read_tsv(Path(a.review) / 's6' / a.tag / ('shard_%03d' % i) / 'occurrences.tsv.gz'):
            if r['occurrence_id'] in need:
                occ[(r['structure'], r['occurrence_id'])] = (r['sample'], r['contig'], r['start0'], r['end0'],
                                                             r['strand'])
    bnd = {}
    for scope in ('same_assembly', 'cross_assembly'):
        for r in uc.read_tsv(Path(a.boundary_dir) / ('boundary_occurrences_%s_%s.tsv.gz' % (scope, a.tag))):
            if r['occurrence_id'] in need:
                bnd[(scope, r['structure'], r['occurrence_id'])] = (r['up_status'], r['down_status'],
                                                                    r['estimated_span'])
    missing = [k for k in (keys or {}) if k not in found]
    rows = []
    for k in list(found) + missing:
        structure, core, scope, oa, ob = k
        f = found.get(k)
        x = dict(pair_id=pair_id(*k), structure=structure, core_key=core, pair_scope=scope)
        x.update(zip(LABELS[4:], f[0] if f else [keys[k].get(c, '') for c in LABELS[4:]]))
        for s, o in (('a', oa), ('b', ob)):
            e = [] if f else ['pair_not_in_boundary_review']
            ro = occ.get((structure, o))
            if ro is None:
                e.append('occurrence_coordinates_missing')
            d = dict(occurrence=o, sample=o.split('|')[1])
            if ro:
                d.update(zip(('sample', 'contig', 'start0', 'end0', 'strand'), ro))
            d['assembly'] = man.get(d['sample'], {}).get('assembly', '')
            if not d['assembly']:
                e.append('sample_not_in_manifest')
            if f and ro:
                src = ('b' if s == 'a' else 'a') if f[1] else s
                up, down = (f[2], f[3]) if src == 'a' else (f[4], f[5])
                d.update(ext_up=up, ext_down=down)
                if up != '' and down != '':
                    st, en = int(ro[2]), int(ro[3])
                    lo, hi = (st - int(up), en + int(down)) if ro[4] == '+' else (st - int(down), en + int(up))
                    d.update(elem_lo=lo, elem_hi=hi)
            b = bnd.get((scope, structure, o), ('', '', ''))
            d.update(boundary_up=b[0], boundary_down=b[1], estimated_span=b[2], input_error=','.join(e))
            x.update({'%s_%s' % (s, q): v for q, v in d.items()})
        rows.append([x.get(c, '') for c in INPUT_COLS])
    return rows


def classify(a, rows, man, info):
    """new / changed / reused per pair against an earlier result's fingerprint index; retired = earlier IDs not in
    this catalogue. Only new and changed pairs are run; reused pairs keep their earlier rows (ucd_combine)."""
    prev = {}
    with gzip.open(a.previous_index, 'rt') as f:
        f.readline()
        for line in f:
            pid, kind, fp, cv = line.rstrip('\n').split('\t')
            prev.setdefault(pid, {})[kind] = (fp, cv)
    cfg, ccfg, code = dict(uc.CONFIG), ck.CHECK_CONFIG, code_hashes()
    keep, cls, n = [], [], Counter()
    for r in rows:
        d = dict(zip(INPUT_COLS, r))
        fps = {k: fingerprint(d, k, man, cfg, ccfg, code) for k in KINDS}
        old = prev.get(r[0])
        if old is None:
            c = 'new'
        elif all(old.get(k, ('',))[0] == fps[k] for k in KINDS):
            c = 'reused'
        else:
            c = 'changed'
        n[c] += 1
        cls.append(dict(pair_id=r[0], status=c, **{'fp_' + k: fps[k] for k in KINDS},
                        **{'previous_fp_' + k: (old or {}).get(k, ('',))[0] for k in KINDS},
                        previous_code_version=(old or {}).get(KINDS[0], ('', ''))[1]))
        if c != 'reused':
            keep.append(r)
    ids = {r[0] for r in rows}
    retired = sorted(set(prev) - ids)
    n['retired'] = len(retired)
    w = Path(a.work)
    uc.write_tsv(w / 'pair_status.tsv.gz', cls + [dict(pair_id=i, status='retired') for i in retired],
                 ['pair_id', 'status'] + ['fp_' + k for k in KINDS] + ['previous_fp_' + k for k in KINDS] +
                 ['previous_code_version'])
    with gzip.open(w / 'reuse_ids.txt.gz', 'wt') as f:
        f.write('pair_id\n' + ''.join(c['pair_id'] + '\n' for c in cls if c['status'] == 'reused'))
    info.update(previous_index=a.previous_index, previous_index_sha256=sha_file(a.previous_index),
                pair_status=dict(n), pair_status_table=str(w / 'pair_status.tsv.gz'),
                pair_status_sha256=sha_file(w / 'pair_status.tsv.gz'), reuse_ids=str(w / 'reuse_ids.txt.gz'),
                decision_version=DECISION_VERSION)
    return keep, prev


def cmd_index(a):
    """Fingerprint index of an earlier final result: every (pair_id, kind) row of <final>/all_rows_<kind>.tsv.gz,
    fingerprinted with the input row and code hashes of the run that produced it (by code_version)."""
    set_kinds(a.kinds.split(','))
    srcs = {}
    for spec in a.source:
        label, inputs, frz = spec.split('|')
        fj = json.loads(Path(frz).read_text())
        srcs[label] = dict(inputs=inputs.split(','), code=fj['code_sha256'], cfg=fj['config'],
                           ccfg=fj.get('check_config', {'max_chain_overlap': 20}))
    man = uc.manifest(a.manifest)
    rows = {}
    for label, sdef in srcs.items():
        for p in sdef['inputs']:
            for r in uc.read_tsv(p):
                rows.setdefault(label, {}).setdefault(r['pair_id'], r)
    out = Path(a.out)
    n = Counter()
    with gzip.open(out, 'wt') as fo:
        fo.write('pair_id\tinterval_kind\tfingerprint\tcode_version\n')
        for kind in KINDS:
            with gzip.open(Path(a.final) / ('all_rows_%s.tsv.gz' % kind), 'rt') as f:
                h = f.readline().rstrip('\n').split('\t')
                ip, iv = h.index('pair_id'), h.index('code_version')
                for line in f:
                    v = line.rstrip('\n').split('\t')
                    pid, cv = v[ip], v[iv]
                    s = srcs[cv]
                    fp = fingerprint(rows[cv][pid], kind, man, s['cfg'], s['ccfg'], s['code'])
                    fo.write('%s\t%s\t%s\t%s\n' % (pid, kind, fp, cv))
                    n[(kind, cv)] += 1
    print(json.dumps({'%s|%s' % k: v for k, v in n.items()}, indent=1))


def cmd_freeze(a):
    out, work = Path(a.out), Path(a.work)
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    man = uc.manifest(a.manifest)
    ci = {c: i for i, c in enumerate(INPUT_COLS)}
    set_kinds(a.kinds.split(','))
    info = {}
    if a.from_input:
        want = ids_of(a.ids)
        rows = [[r[c] for c in INPUT_COLS] for r in uc.read_tsv(a.from_input) if r['pair_id'] in want]
        info.update(from_input=a.from_input, from_input_sha256=sha_file(a.from_input), ids=a.ids,
                    ids_sha256=sha_file(a.ids), ids_requested=len(want), ids_found=len(rows))
        if len(rows) != len(want):
            sys.exit('%d requested IDs not in %s' % (len(want) - len(rows), a.from_input))
    else:
        if a.input_table:  # pre-built pair table in INPUT_COLS (e.g. dna_discovery.py collect)
            rows = [[r.get(c, '') for c in INPUT_COLS] for r in uc.read_tsv(a.input_table)]
            info.update(input_table=a.input_table, input_table_sha256=sha_file(a.input_table))
        else:
            rows = build_from_catalogue(a, man)
        ids = [r[0] for r in rows]
        if len(set(ids)) != len(ids):
            sys.exit('pair_id collision: %d rows, %d distinct IDs' % (len(ids), len(set(ids))))
        info['catalogue_pairs'] = len(rows)
        cat = work / 'catalogue_ids.txt.gz'
        with gzip.open(cat, 'wt') as f:
            f.write('pair_id\n' + '\n'.join(sorted(ids)) + '\n')
        info.update(catalogue_ids=str(cat), catalogue_ids_sha256=sha_file(cat))
        if a.previous_index:
            rows, prev = classify(a, rows, man, info)
        if a.exclude_input:
            ex = ids_of(a.exclude_input)
            rows = [r for r in rows if r[0] not in ex]
            info.update(exclude_input=a.exclude_input, exclude_input_sha256=sha_file(a.exclude_input),
                        excluded_ids=len(ex), excluded_ids_in_catalogue=len(ex & set(ids)),
                        remaining_after_difference=len(rows))
    rows.sort(key=lambda r: (r[ci['a_sample']], r[ci['b_sample']], r[0]))
    n_chunks = (len(rows) + a.chunk_size - 1) // a.chunk_size
    (work / 'inputs').mkdir(exist_ok=True)
    allp = work / 'input_pairs.tsv.gz'
    with gzip.open(allp, 'wt') as fa:
        fa.write('\t'.join(INPUT_COLS) + '\n')
        for c in range(n_chunks):
            part = rows[c * a.chunk_size:(c + 1) * a.chunk_size]
            with gzip.open(work / 'inputs' / ('c%04d.tsv.gz' % c), 'wt') as f:
                f.write('\t'.join(INPUT_COLS) + '\n')
                for r in part:
                    r[ci['chunk']] = c
                    line = '\t'.join(str(v) for v in r) + '\n'
                    f.write(line)
                    fa.write(line)
    samples = sorted({r[ci[s + '_sample']] for r in rows for s in 'ab'})
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
    errs = Counter(e for r in rows for s in 'ab' for e in r[ci[s + '_input_error']].split(',') if e)
    frz = dict(version=VERSION, screen_version=uc.VERSION, check_version=ck.CHECK_VERSION, config=dict(uc.CONFIG),
               check_config=ck.CHECK_CONFIG, coordinates=ck.COORDS, code_sha256=code_hashes(),
               scope=('two checks of a pre-built pair table (%s)' % a.input_table if a.input_table else
                      'full comparison of an existing candidate catalogue (%s, %s); not a complete discovery over all '
                      'FASTA' % ('boundary review core pairs' if a.catalogue == 'boundary' else a.catalogue, a.tag)),
               kinds=list(KINDS),
               inputs={str(p): sha_file(p) for p in [Path(a.manifest)] +
                       ([Path(a.catalogue)] if a.catalogue not in ('boundary', '') else [])},
               input_pairs=len(rows), input_pairs_path=str(allp), input_pairs_sha256=sha_file(allp),
               sources_sha256=sha_file(out / 'sources.tsv.gz'), chunk_size=a.chunk_size, chunks=n_chunks,
               input_errors=dict(errs), sources=dict(Counter(s['status'] for s in src)), work=str(work),
               frozen_at=time.strftime('%Y-%m-%d %H:%M:%S'), **info)
    (out / 'freeze.json').write_text(json.dumps(frz, indent=1) + '\n')
    print(json.dumps({k: v for k, v in frz.items() if k != 'config'}, indent=1))


# ---------------------------------------------------------------- chunk

_W = {}


def _init(man, cfg, ccfg, code_version):
    _W.update(man=man, cfg=cfg, ccfg=ccfg, code_version=code_version, genomes=OrderedDict(), sources=OrderedDict())


def _lru(cache, key, make, size):
    if key in cache:
        cache.move_to_end(key)
        return cache[key]
    v = cache[key] = make()
    if len(cache) > size:
        cache.popitem(last=False)
    return v


def screen_genome(sample):
    return _lru(_W['genomes'], sample, lambda: uc.read_fasta(_W['man'][sample]['fasta']), 64)


def source(sample):
    m = _W['man'][sample]
    return _lru(_W['sources'], sample, lambda: ck.Source(m['fasta'], m['fna_sha256']), 400)


def interval(r, s, kind):
    """(lo, hi), or the reason there is none (input error, or which boundary / extension is missing)."""
    if r[s + '_input_error']:
        return 'input_error:' + r[s + '_input_error']
    if kind == 'cds_core':
        return int(r[s + '_start0']), int(r[s + '_end0'])
    why = []
    for end in ('up', 'down'):
        st = r['%s_boundary_%s' % (s, end)]
        if st != 'resolved':
            why.append('boundary_%s=%s' % (end, st or 'none'))
    if r[s + '_elem_lo'] == '':
        why.append('pair_extension_not_determined')
    if why:
        return 'candidate_interval_not_available(%s)' % ','.join(why)
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
            row = dict({k: r[k] for k in LABELS}, pair_key=key, interval_kind=kind, code_version=_W['code_version'])
            for s in 'ab':
                row.update({'%s_%s' % (s, f): r[s + '_' + f] for f in ('occurrence', 'sample', 'assembly', 'contig',
                                                                         'strand', 'boundary_up', 'boundary_down',
                                                                         'estimated_span')})
            ivs = {s: interval(r, s, kind) for s in 'ab'}
            bad = [s for s in 'ab' if not isinstance(ivs[s], tuple)]
            if bad:
                why = ';'.join('%s:%s' % (s, ivs[s]) for s in bad)
                row.update(source_check='not_run', source_reasons=why, ucd_check='not_run', final_status='not_passed',
                           failure_reason=why)
                rows.append(row)
                continue
            # screen (its own extraction)
            S = {}
            for s in 'ab':
                cs = screen_genome(r[s + '_sample']).get(r[s + '_contig'])
                S[s] = None if cs is None else dict(uc.sides_of(cs, ivs[s][0], ivs[s][1], r[s + '_strand'], L),
                                                    id=r[s + '_occurrence'], sample=r[s + '_sample'],
                                                    contig=r[s + '_contig'])
            if S['a'] and S['b']:
                sc = uc.compare(S['a'], S['b'], cfg, tmp)
                row.update({'screen_' + k: sc[k] for k in SCREEN})
            else:
                row['screen_verdict'] = 'screen_error:contig_missing'
            # Check 1
            C, why = {}, []
            for s in 'ab':
                src = source(r[s + '_sample'])
                c = C[s] = ck.check1_side(src, r[s + '_contig'], ivs[s][0], ivs[s][1], r[s + '_strand'], L, S[s])
                row.update({'%s_fasta_sha256' % s: getattr(src, 'sha256', ''), '%s_contig_len' % s: c['contig_len'],
                            '%s_edge' % s: c['edge']})
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
    set_kinds(frz.get('kinds', uc.KINDS))
    if code_hashes() != frz['code_sha256']:
        sys.exit('code differs from the frozen code (freeze.json); refusing to run')
    idx = a.index if a.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID']) + a.offset
    if idx >= frz['chunks']:
        print('chunk %d beyond %d chunks; nothing to do' % (idx, frz['chunks']))
        return
    d = work / 'chunks' / ('c%04d' % idx)
    done = d / 'done.json'
    if done.exists():
        j = json.loads(done.read_text())
        if j.get('input_pairs_sha256') == frz['input_pairs_sha256'] and \
                sha_file(d / 'rows.tsv.gz') == j.get('rows_sha256'):
            print('chunk %d already complete; skipped' % idx)
            return
    pairs = uc.read_tsv(work / 'inputs' / ('c%04d.tsv.gz' % idx))
    man = uc.manifest(a.manifest)
    cv = '%s+%s' % (VERSION, hashlib.sha256(''.join(frz['code_sha256'][c] for c in CODE).encode()).hexdigest()[:12])
    d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows, seqs, hsps, cpu, wall = [], {}, [], [], []
    with ProcessPoolExecutor(a.threads, initializer=_init,
                             initargs=(man, frz['config'], frz['check_config'], cv)) as ex:
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
        rows_sha256=sha_file(d / 'rows.tsv.gz'), input_pairs_sha256=frz['input_pairs_sha256'], code_version=cv,
        wall_s=round(time.time() - t0, 1), threads=a.threads, worker_wall_s=round(sum(wall), 1),
        child_process_cpu_s=round(sum(cpu), 1), host=socket.gethostname(),
        slurm=dict(job=os.environ.get('SLURM_ARRAY_JOB_ID', os.environ.get('SLURM_JOB_ID', '')),
                   task=os.environ.get('SLURM_ARRAY_TASK_ID', ''), restart=os.environ.get('SLURM_RESTART_COUNT', '')),
        finished=time.strftime('%Y-%m-%d %H:%M:%S'))
    done.write_text(json.dumps(j, indent=1) + '\n')
    print(json.dumps(j))


# ---------------------------------------------------------------- merge (one run; streamed)

def cmd_merge(a):
    out, work = Path(a.out), Path(a.work)
    frz = json.loads((out / 'freeze.json').read_text())
    set_kinds(frz.get('kinds', uc.KINDS))
    expect = set()
    for i in range(frz['chunks']):
        expect |= {r['pair_id'] for r in uc.read_tsv(work / 'inputs' / ('c%04d.tsv.gz' % i))}
    missing_chunks, done = [], []
    for i in range(frz['chunks']):
        d = work / 'chunks' / ('c%04d' % i)
        ok = (d / 'done.json').exists()
        if ok:
            j = json.loads((d / 'done.json').read_text())
            ok = sha_file(d / 'rows.tsv.gz') == j['rows_sha256']
        if ok:
            done.append(j)
        else:
            missing_chunks.append(i)
    if missing_chunks:
        comp = dict(complete=False, missing_chunks=missing_chunks)
        (out / 'completeness.json').write_text(json.dumps(comp, indent=1) + '\n')
        sys.exit('incomplete: %d chunks missing; resubmit them (finished chunks are skipped)' % len(missing_chunks))
    seen = {k: set() for k in KINDS}
    dup, extra = Counter(), Counter()
    outs = {k: gzip.open(work / ('rows_%s.tsv.gz' % k), 'wt') for k in KINDS}
    for f in outs.values():
        f.write('\t'.join(ROW_COLS) + '\n')
    for i in range(frz['chunks']):
        with gzip.open(work / 'chunks' / ('c%04d' % i) / 'rows.tsv.gz', 'rt') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            assert hdr == ROW_COLS
            ip, ik = hdr.index('pair_id'), hdr.index('interval_kind')
            for line in f:
                v = line.split('\t')
                pid, kind = v[ip], v[ik]
                if pid not in expect or kind not in seen:
                    extra[kind] += 1
                elif pid in seen[kind]:
                    dup[kind] += 1
                seen.setdefault(kind, set()).add(pid)
                outs[kind].write(line) if kind in outs else None
    for f in outs.values():
        f.close()
    missing = {k: len(expect - seen[k]) for k in KINDS}
    comp = dict(expected_pairs=len(expect), expected_rows=len(expect) * len(KINDS),
                observed_rows={k: len(seen[k]) for k in KINDS}, missing=missing, duplicated=dict(dup),
                unexpected=dict(extra), missing_chunks=[],
                complete=not (any(missing.values()) or dup or extra))
    (out / 'completeness.json').write_text(json.dumps(comp, indent=1) + '\n')
    if not comp['complete']:
        sys.exit('incomplete or inconsistent: %s' % json.dumps(comp))
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
        fo.write('\t'.join(HSP_COLS) + '\n')
        for i in range(frz['chunks']):
            with gzip.open(work / 'chunks' / ('c%04d' % i) / 'hsps.tsv.gz', 'rt') as f:
                f.readline()
                shutil.copyfileobj(f, fo)
    (out / 'chunks_done.json').write_text(json.dumps(done, indent=1) + '\n')
    summ = dict(version=VERSION, completeness=comp, code_versions=dict(Counter(j['code_version'] for j in done)),
                timing=dict(chunks=len(done), sum_chunk_wall_s=round(sum(j['wall_s'] for j in done), 1),
                            worker_cpu_s=round(sum(j['child_process_cpu_s'] for j in done), 1),
                            worker_wall_s=round(sum(j['worker_wall_s'] for j in done), 1),
                            first_finished=min(j['finished'] for j in done),
                            last_finished=max(j['finished'] for j in done)),
                evidence_sequences_distinct=len(seen_s),
                files={p.name: dict(path=str(p), bytes=p.stat().st_size, sha256=sha_file(p)) for p in
                       [work / ('rows_%s.tsv.gz' % k) for k in KINDS] + [seq_p, hsp_p]})
    (out / 'merge.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


# ---------------------------------------------------------------- cleanup

def cmd_cleanup(a):
    out, work = Path(a.out), Path(a.work).resolve()
    comp = json.loads((out / 'completeness.json').read_text())
    summ = json.loads((out / 'merge.json').read_text())
    if not comp.get('complete'):
        sys.exit('merge not complete; nothing removed')
    for name, f in summ['files'].items():
        if sha_file(f['path']) != f['sha256']:
            sys.exit('merged file %s does not match merge.json; nothing removed' % name)
    targets = [work / 'chunks', work / 'inputs']
    for t in targets:
        if not str(t).startswith('/global/scratch/users/kh36969/longread_callout/core_context/') or t.parent != work:
            sys.exit('refusing to remove %s (outside this run directory)' % t)
    listing = [dict(path=str(p), bytes=p.stat().st_size) for t in targets if t.exists() for p in sorted(t.rglob('*'))
               if p.is_file()]
    freed = sum(x['bytes'] for x in listing)
    uc.write_tsv(out / 'cleanup_list.tsv.gz', listing, ['path', 'bytes'])
    if a.apply:
        for t in targets:
            if t.exists():
                shutil.rmtree(t)
    rec = dict(removed=[str(t) for t in targets], files=len(listing), bytes=freed, applied=a.apply,
               list='cleanup_list.tsv.gz', kept=sorted(summ['files']) + ['input_pairs.tsv.gz', 'sources.tsv.gz',
                                                                         'freeze.json'],
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
            x.add_argument('--catalogue', default='',
                           help="'boundary' (every boundary-review core pair) or a catalogue file")
            x.add_argument('--exclude-input', default='', help='drop the pair IDs of this earlier input table')
            x.add_argument('--from-input', default='', help='take rows from this earlier frozen input table ...')
            x.add_argument('--ids', default='', help='... for these pair IDs (one per line or a table with pair_id)')
            x.add_argument('--review', default=str(uc.SCRATCH / 'scale_comamonadaceae_v0.1' / 'review'))
            x.add_argument('--tag', default='b1-2-3')
            x.add_argument('--boundary-dir', default=str(uc.REPO), help='dir with boundary_occurrences_<scope>_<tag>')
            x.add_argument('--previous-index', default='', help='fingerprint index of an earlier result (index)')
            x.add_argument('--input-table', default='', help='pre-built input table (INPUT_COLS); all its rows')
            x.add_argument('--kinds', default=','.join(uc.KINDS), help='interval kinds of this run')
            x.add_argument('--chunk-size', type=int, default=5000)
        if name == 'chunk':
            x.add_argument('--index', type=int)
            x.add_argument('--offset', type=int, default=0, help='added to SLURM_ARRAY_TASK_ID (arrays > 1000)')
            x.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
        if name == 'cleanup':
            x.add_argument('--apply', action='store_true')
    x = sp.add_parser('index')
    x.add_argument('--final', required=True, help='dir with all_rows_<kind>.tsv.gz')
    x.add_argument('--kinds', default=','.join(uc.KINDS), help='interval kinds of the final result')
    x.add_argument('--source', action='append', required=True,
                   help="'<code_version>|<input_pairs[,input_pairs]>|<freeze.json>' of a run whose rows are in --final")
    x.add_argument('--manifest', default=str(uc.REPO / 'manifest.tsv'))
    x.add_argument('--out', required=True)
    a = p.parse_args()
    if a.cmd == 'freeze' and not (a.catalogue or a.from_input or a.input_table):
        sys.exit('freeze needs --catalogue or --from-input/--ids or --input-table')
    dict(freeze=cmd_freeze, chunk=cmd_chunk, merge=cmd_merge, cleanup=cmd_cleanup, index=cmd_index)[a.cmd](a)


if __name__ == '__main__':
    main()
