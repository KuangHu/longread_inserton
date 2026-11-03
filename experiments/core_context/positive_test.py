"""Computational positive test of the read-site search (fixed before running; reported apart from real samples).

prepare  In validation samples (never development samples), seeded random catalog occurrences of panel cores are cut
         out of the sample's own assembly together with an extra 0-300 bp of host sequence per side. The real reads
         still carry the element, so at each cut the search should report a position not in the (edited) assembly,
         whose insertion = extra left + core + extra right. Strata (each cell drawn apart, seed SEED):
           core length   < 1000 / 1000-1499 / >= 1500 bp (coding span of the occurrence)
           structure     single_gene / two_gene
           repeat        the core's copies in the sample's catalog: 1 (single copy) / >= 2 (repeated background)
           extra         none (0 + 0 bp) / extra (each side drawn uniformly from 1-300 bp)
         Allowed cuts: extra + core + extra < MAX_INSERT; no other catalog occurrence of any core lies partly in the
         cut window (occurrences wholly inside, e.g. the genes of a two-gene core, go with the cut; the catalog has
         one at almost every gene, so no margin there); cuts are >= MARGIN apart and >= MARGIN
         from a contig end. At most PER_SAMPLE
         cuts per sample; up to PER_CELL per cell.
         Inputs of the edited run (all in edited coordinates, under <work>/edited, read by read_site_search via
         set_work): edited assemblies, manifest.tsv (fasta, total_bp, sha256), genes/<s>.tsv.gz (genes overlapping a
         cut dropped, others shifted, contig_length updated) and occurrences.tsv (the sample's catalog occurrences, cut
         ones dropped, others shifted). The panel is built from the edited catalog, so the cut copy itself is not a
         panel sequence. The truth (original and edited coordinates, extra lengths, cut sequence) is written only
         to --truth, outside <work>; no search stage reads it.
score    Per cut, the merged positions of the edited run at the main setting whose two junctions both lie within
         TOL of the cut junction: detected (not_in_assembly with >= MAIN_COUNT molecules), detected_below_threshold,
         other_class, or not_detected; insertion length and extra per side against the truth; boundary class; size
         class. Whether the same occurrence was found as a known position in the real (unedited) search of the
         sample is reported beside it (separates read coverage from the method). Not-in-assembly positions of the
         edited run away from every cut are compared with the real run of the same sample (new ones = false
         positives of the edit test). With the candidate review (structure_summary.tsv), the consensus insertion
         is compared with the cut sequence.
readmap  One edited sample (--samples s): its reads against its edited assembly (the genome-stage alignment the
         empty-site count and the candidate evidence read), under <work>/edited/readmap/<s>.
Usage: python positive_test.py prepare --samples a,b --work <scratch> --truth <dir>
       python positive_test.py score --out <edited batch repo dir> --work <scratch> --truth <dir> --real <dir,...>
"""
import argparse
import gzip
import hashlib
import json
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import protein_context as pc
import read_site_search as rs
import s3_readcheck as rc3
import ucd_readcheck as u

SEED = 20261102
EXTRA_MAX = 300
MARGIN = 10000  # bp of untouched sequence around a cut window, and from contig ends
PER_SAMPLE = 10
PER_CELL = 6
LEN_BINS = ((0, 1000, '<1000'), (1000, 1500, '1000-1499'), (1500, 10 ** 9, '>=1500'))
OCC = u.PROJ / 'core_context' / 'scale_comamonadaceae_v0.1' / 'review' / 's6' / 'b1-2-3'


def len_bin(n):
    return next(name for lo, hi, name in LEN_BINS if lo <= n < hi)


def occurrences(samples):
    out = defaultdict(list)
    for k in range(64):
        for r in pc.read_tsv(OCC / ('shard_%03d' % k) / 'occurrences.tsv.gz'):
            if r['sample'] in samples:
                out[r['sample']].append(r)
    return out


def panel_cores():
    return {(p['structure'], p['core_key']) for p in pc.read_tsv(rs.ROUND5 / 'pairs.tsv.gz')
            if p['tier'] in ('A', 'B') and u.pair_role(p) == 'candidate' and p['pair_class'] == 'both_positions_supported'}


def shift(pos, cuts):
    """Edited coordinate of an original position outside every cut window (cuts: sorted (lo, hi) of one contig)."""
    return pos - sum(hi - lo for lo, hi in cuts if hi <= pos)


def inside(lo, hi, cuts):
    return any(lo < b and hi > a for a, b in cuts)


def cmd_prepare(a):
    samples = a.samples.split(',')
    dev = set()
    for d in ('read_site_search_pilot',):
        dev |= set(rs.pilot_samples(rs.REPO / d))
    assert not dev & set(samples), ('development samples', dev & set(samples))
    man = rs.manifest()
    cores = panel_cores()
    occ = occurrences(set(samples))
    rng = random.Random(SEED)
    seqs = {s: rc3.read_fasta(man[s]['fasta']) for s in samples}
    cand = defaultdict(list)  # cell -> occurrences
    for s in samples:
        n = Counter(r['core_key'] for r in occ[s])
        for r in occ[s]:
            if (r['structure'], r['core_key']) not in cores:
                continue
            L = int(r['end0']) - int(r['start0'])
            for extra in ('none', 'extra'):
                cell = (len_bin(L), r['structure'], 'single_copy' if n[r['core_key']] == 1 else 'repeated', extra)
                cand[cell].append(r)
    cuts, used = defaultdict(list), Counter()  # sample -> chosen
    chosen_occ = set()
    for cell in sorted(cand):
        pool = sorted(cand[cell], key=lambda r: r['occurrence_id'])
        rng.shuffle(pool)
        got = 0
        for r in pool:
            if got >= PER_CELL:
                break
            s = r['sample']
            if used[s] >= PER_SAMPLE or r['occurrence_id'] in chosen_occ:
                continue
            el, er = (0, 0) if cell[3] == 'none' else (rng.randint(1, EXTRA_MAX), rng.randint(1, EXTRA_MAX))
            lo, hi = int(r['start0']) - el, int(r['end0']) + er
            clen = len(seqs[s][r['contig']])
            if hi - lo >= rs.MAX_INSERT or lo < MARGIN or hi > clen - MARGIN:
                continue
            if any(o['contig'] == r['contig'] and o['occurrence_id'] != r['occurrence_id']
                   and int(o['start0']) < hi and int(o['end0']) > lo
                   and not (int(o['start0']) >= lo and int(o['end0']) <= hi) for o in occ[s]):
                continue
            if any(c['contig'] == r['contig'] and c['lo'] < hi + MARGIN and c['hi'] > lo - MARGIN for c in cuts[s]):
                continue
            cuts[s].append(dict(occurrence_id=r['occurrence_id'], sample=s, contig=r['contig'], structure=r['structure'],
                                core_key=r['core_key'], strand=r['strand'], core_start0=int(r['start0']),
                                core_end0=int(r['end0']), core_len=int(r['end0']) - int(r['start0']), extra_left=el,
                                extra_right=er, lo=lo, hi=hi, insert_len=hi - lo, len_bin=cell[0],
                                repeat=cell[2], extra_mode=cell[3]))
            chosen_occ.add(r['occurrence_id'])
            used[s] += 1
            got += 1
    ed = a.work / 'edited'
    (ed / 'genes').mkdir(parents=True, exist_ok=True)
    a.truth.mkdir(parents=True, exist_ok=True)
    assert a.work.resolve() not in (a.truth.resolve(),) + tuple(a.truth.resolve().parents), 'truth must be outside the work directory'
    mrows, orows, truth = [], [], []
    for s in samples:
        by = defaultdict(list)
        for c in cuts[s]:
            by[c['contig']].append((c['lo'], c['hi']))
        for v in by.values():
            v.sort()
        fa = ed / ('%s.edited.fna' % s)
        with open(str(fa), 'w') as h:
            for name, sq in seqs[s].items():
                parts, last = [], 0
                for lo, hi in by.get(name, []):
                    parts.append(sq[last:lo])
                    last = hi
                parts.append(sq[last:])
                e = ''.join(parts)
                h.write('>%s\n' % name + '\n'.join(e[i:i + 80] for i in range(0, len(e), 80)) + '\n')
        total = sum(len(x) for x in rc3.read_fasta(str(fa)).values())
        mrows.append(dict(sample=s, fasta=str(fa), total_bp=total,
                          fna_sha256=hashlib.sha256(fa.read_bytes()).hexdigest(), edited_from=man[s]['fasta'],
                          cuts=len(cuts[s])))
        for r in occ[s]:
            lo, hi = int(r['start0']), int(r['end0'])
            if inside(lo, hi, by.get(r['contig'], [])):
                continue
            cs = by.get(r['contig'], [])
            orows.append(dict(r, start0=shift(lo, cs), end0=shift(hi, cs)))
        grows = []
        for g in pc.read_tsv(rs.GENES / ('%s.tsv.gz' % s)):
            cs = by.get(g['contig'], [])
            lo, hi = int(g['start0']), int(g['end0'])
            if inside(lo, hi, cs):
                continue
            grows.append(dict(g, start0=shift(lo, cs), end0=shift(hi, cs),
                              contig_length=int(g['contig_length']) - sum(y - x for x, y in cs)))
        if grows:
            pc.write_tsv(ed / 'genes' / ('%s.tsv.gz' % s), grows, list(grows[0]))
        for c in cuts[s]:
            cs = by[c['contig']]
            cut_seq = seqs[s][c['contig']][c['lo']:c['hi']].upper()
            truth.append(dict(c, edited_junction=shift(c['lo'], cs),
                              cut_sha256=hashlib.sha256(cut_seq.encode()).hexdigest()))
            with open(str(a.truth / 'cut_sequences.fa'), 'a') as h:
                h.write('>%s|%s:%d-%d\n%s\n' % (c['occurrence_id'], c['contig'], c['lo'], c['hi'], cut_seq))
    pc.write_tsv(ed / 'manifest.tsv', mrows, list(mrows[0]))
    pc.write_tsv(ed / 'occurrences.tsv', orows, list(orows[0]))
    pc.write_tsv(a.truth / 'truth.tsv', truth, list(truth[0]))
    summ = dict(seed=SEED, samples=samples, cuts=len(truth), per_sample=dict(Counter(t['sample'] for t in truth)),
                cells={'|'.join(k): v for k, v in sorted(Counter((t['len_bin'], t['structure'], t['repeat'],
                                                                   t['extra_mode']) for t in truth).items())},
                cells_available={'|'.join(k): len(v) for k, v in sorted(cand.items())},
                per_sample_max=PER_SAMPLE, per_cell_max=PER_CELL, extra_max=EXTRA_MAX, margin=MARGIN,
                edited_manifest_sha256=hashlib.sha256((ed / 'manifest.tsv').read_bytes()).hexdigest(),
                truth_sha256=hashlib.sha256((a.truth / 'truth.tsv').read_bytes()).hexdigest())
    (a.truth / 'prepare.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


def cmd_score(a):
    truth = list(pc.read_tsv(a.truth / 'truth.tsv'))
    rows = [r for r in pc.read_tsv(a.out / 'positions.tsv.gz') if int(r['flank']) == rs.MAIN_FLANK and r['kind'] == 'core']
    real = []
    for d in a.real.split(','):
        real += [r for r in pc.read_tsv(Path(d) / 'positions.tsv.gz') if int(r['flank']) == rs.MAIN_FLANK
                 and r['kind'] == 'core']
    found_real = {o for r in real if r['cls'].startswith('known_catalog_position') for o in r['catalog'].split(';')}
    groups = rs.family_groups(a.work)
    st = {}
    sp = a.out / 'candidates' / 'structure_summary.tsv'
    if sp.exists():
        for r in pc.read_tsv(sp):
            st[(r['sample'], r['candidate'].split(':')[1], r['candidate'].split(':')[2])] = r
    out, near_cut = [], set()
    for t in truth:
        J = int(t['edited_junction'])
        hits = [r for r in rows if r['sample'] == t['sample'] and r['left_contig'] == r['right_contig'] == t['contig']
                and r['left_pos'] and abs(int(r['left_pos']) - J) <= rs.TOL and abs(int(r['right_pos']) - J) <= rs.TOL]
        near_cut |= {id(r) for r in hits}
        rel = [(r, rs.core_relation(r['core_key'], t['core_key'], groups)[0]) for r in hits]
        rel = [(r, k) for r, k in rel if k != 'unrelated']
        nia = [r for r, _ in rel if r['cls'].startswith('not_in_assembly')]
        best = max(nia or [r for r, _ in rel], key=lambda r: int(r['molecules']), default=None)
        if best is None:
            res = 'not_detected'
        elif best['cls'].startswith('not_in_assembly'):
            res = 'detected' if int(best['molecules']) >= rs.MAIN_COUNT else 'detected_below_threshold'
        else:
            res = 'other_class:' + best['cls']
        e = dict(t, result=res, real_run_known=t['occurrence_id'] in found_real, seeds=len(rel))
        if best is not None:
            e.update(cls=best['cls'], molecules=best['molecules'], core_key_found=best['core_key'],
                     boundary=best['boundary'], insert_len_found=best['insert_len'],
                     insert_len_error=int(best['insert_len']) - int(t['insert_len']) if best['insert_len'] else '',
                     extra_left_found=best['extra_left'], extra_right_found=best['extra_right'],
                     extra_error=max(abs(int(float(best['extra_left'])) - int(t['extra_left'])),
                                     abs(int(float(best['extra_right'])) - int(t['extra_right'])))
                     if best['extra_left'] != '' else '',
                     size_class_found=best['insert_size'], size_class_true=rs.size_class(int(t['insert_len'])))
        s_ = [v for (s, c, jj), v in st.items() if s == t['sample'] and c == t['contig']
              and abs(int(jj.split('-')[0]) - J) <= rs.TOL]
        if s_:
            v = s_[0]
            e.update(consensus_status=v['status'], consensus_insert_len=v.get('insert_len', ''),
                     consensus_target=v.get('target_1_2_genes', ''), consensus_genes=v.get('gene_count_range', ''))
            if v.get('insert_len'):
                e['consensus_len_error'] = int(v['insert_len']) - int(t['insert_len'])
        out.append(e)
    pc.write_tsv(a.out / 'positive_test_scores.tsv', out, sorted({k for e in out for k in e},
                                                                    key=lambda k: list(out[0]).index(k)
                                                                    if k in out[0] else 999))
    realnia = {(r['sample'], r['left_contig'], int(r['left_pos']) // rs.TOL) for r in real
               if r['cls'].startswith('not_in_assembly') and r['left_pos']}
    away = [r for r in rows if id(r) not in near_cut and r['cls'].startswith('not_in_assembly')
            and int(r['molecules']) >= rs.MAIN_COUNT]

    def strata(key):
        c = defaultdict(Counter)
        for e in out:
            c[e[key]][e['result'].split(':')[0]] += 1
        return {k: dict(v) for k, v in sorted(c.items())}
    det = [e for e in out if e['result'] == 'detected']
    summ = dict(
        cuts=len(out), result=dict(Counter(e['result'].split(':')[0] for e in out)),
        other_classes=dict(Counter(e['result'] for e in out if e['result'].startswith('other_class'))),
        among_real_run_known=dict(Counter(e['result'].split(':')[0] for e in out if e['real_run_known'])),
        among_not_real_run_known=dict(Counter(e['result'].split(':')[0] for e in out if not e['real_run_known'])),
        by_len_bin=strata('len_bin'), by_structure=strata('structure'), by_repeat=strata('repeat'),
        by_extra_mode=strata('extra_mode'),
        detected_boundary=dict(Counter(e['boundary'] for e in det)),
        detected_insert_len_abs_error=dict(Counter(
            '0-10' if abs(e['insert_len_error']) <= 10 else '11-50' if abs(e['insert_len_error']) <= 50 else '>50'
            for e in det if e.get('insert_len_error') != '')),
        detected_extra_max_side_error=dict(Counter(
            '0-20' if e['extra_error'] <= 20 else '21-50' if e['extra_error'] <= 50 else '>50'
            for e in det if e.get('extra_error') != '')),
        detected_size_class_agrees=sum(e['size_class_found'] == e['size_class_true'] for e in det),
        consensus=dict(Counter('%s|%s' % (e.get('consensus_status', 'not_reviewed'), e.get('consensus_target', ''))
                               for e in det)),
        consensus_len_abs_error=dict(Counter(
            '0-10' if abs(e['consensus_len_error']) <= 10 else '11-50' if abs(e['consensus_len_error']) <= 50 else '>50'
            for e in det if e.get('consensus_len_error', '') != '')),
        not_in_assembly_away_from_cuts=dict(
            total=len(away), also_in_real_run=sum((r['sample'], r['left_contig'], int(r['left_pos']) // rs.TOL)
                                                 in realnia for r in away if r['left_pos'])))
    (a.out / 'positive_test_score.json').write_text(json.dumps(summ, indent=1) + '\n')
    print(json.dumps(summ, indent=1))


def cmd_readmap(a):
    """Every read of every verified run of one edited sample against its edited assembly, with the genome-stage
    minimap2 arguments (ucd_readcheck cmd_genome); written to <work>/edited/readmap/<sample>."""
    import time
    rs.set_work(a.work)
    s = a.samples
    fa = rs.manifest()[s]['fasta']
    od = rs.readmap_dir(s)
    od.mkdir(parents=True, exist_ok=True)
    runs = u.run_info()[s]
    u.reads_ready.require([r['run'] for r in runs])
    done = []
    for r in runs:
        paf = od / ('%s.paf.gz' % r['run'])
        t0 = time.time()
        cmd = 'minimap2 -x %s %s -t %d %s %s | gzip -1 > %s.tmp' % (
            rc3.PRESET[r['read_type']], ' '.join(rc3.MM_GENOME), a.threads, fa, r['local_file'], paf)
        subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], check=True, stderr=subprocess.DEVNULL)
        Path(str(paf) + '.tmp').rename(paf)
        done.append(dict(run=r['run'], status='mapped', preset=rc3.PRESET[r['read_type']], minimap2_args=rc3.MM_GENOME,
                         assembly=fa, sha256=r['sha256'], seconds=round(time.time() - t0)))
    (od / 'genome.json').write_text(json.dumps(dict(sample=s, runs=done), indent=1) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mode', choices=['prepare', 'readmap', 'score'])
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--samples', default='')
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--truth', type=Path)
    p.add_argument('--out', type=Path)
    p.add_argument('--real', default='')
    a = p.parse_args()
    rs.set_work(a.work) if a.mode == 'score' else None
    dict(prepare=cmd_prepare, readmap=cmd_readmap, score=cmd_score)[a.mode](a)


if __name__ == '__main__':
    main()
