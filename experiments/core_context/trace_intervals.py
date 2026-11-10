"""Trace the s3_readcheck v0.3 rules step by step for the intervals of one assembly in one check directory.

Reads the files the check stage left (alleles.fa, alleles_<run>.paf.gz, genome_<run>.paf.gz) and applies the same
functions the caller uses (test_geometry, cigar_ok, read_calls); nothing is re-aligned. Per interval and read, the
step where the read leaves the count is recorded:
  no_allele_hit        near the interval on the genome, no alignment to this interval's P or A allele
  partial:short_read   no allele alignment spans [lo - ANCHOR, hi + ANCHOR] and the read is shorter than that span
  partial:clipped      no spanning alignment although the read is long enough: the alignment stops inside the span
                       and the read continues (clip_bp = read bases beyond the stop)
  cigar_fail           spans, but the CIGAR limit fails on every spanning alignment
  ambiguous            spans both alleles with an AS difference < MIN_DELTA_AS
  off_locus / locus_ambiguous   genome-wide anchor competition
  P / A                counted call; *_matched when the read is at least as long as the tested span
Usage: python trace_intervals.py <check_dir> <positions.tsv> <sample> <out.tsv> [test_id,...]
"""
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

import scale_entry as se  # noqa: F401  (puts fna_pool on sys.path)
import s3_readcheck as rc3
import protein_context as pc
import ucd_readcheck as u


def paf(path):
    with gzip.open(str(path), 'rt') as h:
        for line in h:
            f = line.rstrip('\n').split('\t')
            tags = dict(x.split(':', 2)[::2] for x in f[12:])
            yield dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], tlen=int(f[6]),
                       ts=int(f[7]), te=int(f[8]), score=int(tags.get('AS', 0)), tp=tags.get('tp', 'P'),
                       cg=tags.get('cg', ''))


def main(cd, positions_tsv, sample, out, only=''):
    cd = Path(cd)
    pos = [p for p in pc.read_tsv(positions_tsv) if p['sample'] == sample]
    want = set(only.split(',')) if only else {p['test_id'] for p in pos}
    seqs = rc3.read_fasta(u.fasta_of(Path(u.__file__).resolve().parent / 'scale_comamonadaceae_v0.1', sample))
    tests = {t['test_id']: t for t in u.make_tests(pos, seqs)[0]}
    geo = {k: rc3.test_geometry(t) for k, t in tests.items()}
    rows = []
    for apaf in sorted(cd.glob('alleles_*.paf.gz')):
        run = apaf.name[len('alleles_'):-len('.paf.gz')]
        galns = defaultdict(list)
        for g in paf(cd / ('genome_%s.paf.gz' % run)):
            galns[g['q']].append(g)
        hits = defaultdict(lambda: defaultdict(list))  # read -> (tid, allele) -> alignments
        for b in paf(apaf):
            tid, role, al = b['t'].split('|')
            if tid in want:
                hits[b['q']][(tid, al)].append(b)
        for p in pos:
            tid = p['test_id']
            if tid not in want:
                continue
            lo, hi = int(p['lo']), int(p['hi'])
            g = geo[tid]['C']
            near = {q for q, gs in galns.items()
                    if any(x['t'] == p['contig'] and x['ts'] < hi + rc3.FLANK and x['te'] > lo - rc3.FLANK for x in gs)}
            for q in sorted(near | {q for q in hits if any(k[0] == tid for k in hits[q])}):
                qlen = next(iter(galns[q] or [x for v in hits[q].values() for x in v]))['qlen']
                best, detail = {}, {}
                for al in 'PA':
                    ilo, ihi, _ = g['iv'][al]
                    for b in hits[q].get((tid, al), []):
                        st = (2 if rc3.cigar_ok(b, ilo, ihi, g['limit']) else 1) if b['ts'] <= ilo and b['te'] >= ihi else 0
                        cur = best.get((tid, 'C', al))
                        if cur is None or (st, b['score']) > cur[:2]:
                            best[(tid, 'C', al)] = (st, b['score'], (b['strand'], rc3.q_at(b, ilo), rc3.q_at(b, ihi - 1))
                                                     if st == 2 else None)
                            detail[al] = b
                if not best:
                    step, why = 'no_allele_hit', ''
                else:
                    call = rc3.read_calls(best, galns[q], geo, qlen).get((tid, 'C'))
                    if call:
                        step, why = call, ''
                    else:  # partial on both alleles: is the read too short, or does the alignment stop?
                        al = max(detail, key=lambda a: detail[a]['te'] - detail[a]['ts'])
                        b, (ilo, ihi, _) = detail[al], g['iv'][al]
                        miss_l, miss_r = max(0, b['ts'] - ilo), max(0, ihi - b['te'])
                        left_q = b['qs'] if b['strand'] == '+' else b['qlen'] - b['qe']  # read bases before the start
                        right_q = b['qlen'] - b['qe'] if b['strand'] == '+' else b['qs']
                        clipped = (miss_l and left_q > 100) or (miss_r and right_q > 100)
                        step = 'partial:clipped' if clipped else 'partial:short_read' if qlen < ihi - ilo else 'partial:ends'
                        why = '%s ts=%d te=%d need=%d-%d missing_left=%d missing_right=%d clip_left=%d clip_right=%d' % (
                            al, b['ts'], b['te'], ilo, ihi, miss_l, miss_r, left_q, right_q)
                gen = ';'.join('%s:%d-%d:%s:%s' % (x['t'], x['ts'], x['te'], x['strand'], x['tp'])
                               for x in sorted(galns[q], key=lambda x: x['qs']))
                rows.append(dict(test_id=tid, run=run, read=q, qlen=qlen, step=step, detail=why, genome_alns=gen))
    pc.write_tsv(out, rows, ['test_id', 'run', 'read', 'qlen', 'step', 'detail', 'genome_alns'])
    by = defaultdict(Counter)
    for r in rows:
        by[r['test_id']][r['step']] += 1
    for p in pos:
        if p['test_id'] in want:
            print(p['test_id'], p['contig'], p['lo'], p['hi'], dict(by[p['test_id']]))


if __name__ == '__main__':
    main(*sys.argv[1:])
