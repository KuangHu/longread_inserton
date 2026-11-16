"""Target-length table of a combined ucd result (result filter only; the two checks are not rerun).

Main table: rows that pass both checks (source_check pass and ucd_check ucd_supported) and whose candidate_insertion
X, as actually extracted in Check 1 (x_len_a / x_len_b), has TARGET_MIN <= length < TARGET_MAX bp on BOTH sides.
The wider strict table (0 < length < 2500) is kept unchanged for traceability; this table is a subset of it.
Counts are pairs, distinct X intervals (sample, contig, interval), assemblies, and same / cross assembly, kept apart.
"""
import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

TARGET_MIN, TARGET_MAX = 700, 2500
KIND = 'candidate_insertion'


def in_target(n):
    return TARGET_MIN <= n < TARGET_MAX


def summarize(rows):
    iv, asm, scope = set(), set(), Counter()
    for r in rows:
        scope[r['pair_scope']] += 1
        for s in 'ab':
            iv.add((r[s + '_sample'], r[s + '_contig'], r[s + '_x_interval']))
            asm.add(r[s + '_sample'])
    return dict(pairs=sum(scope.values()), distinct_x_intervals=len(iv), assemblies=len(asm), by_scope=dict(scope))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--final', required=True, help='dir with pass_/strict_lt2500_candidate_insertion.tsv.gz')
    ap.add_argument('--out', default='', help='output dir (default: --final)')
    a = ap.parse_args()
    fin = Path(a.final)
    out = Path(a.out or a.final)
    out.mkdir(parents=True, exist_ok=True)
    name = 'strict_%d_%d_%s.tsv.gz' % (TARGET_MIN, TARGET_MAX, KIND)
    wide, tgt, below = [], [], Counter()
    with gzip.open(fin / ('strict_lt2500_%s.tsv.gz' % KIND), 'rt') as f, gzip.open(out / name, 'wt') as fo:
        hdr = f.readline()
        h = hdr.rstrip('\n').split('\t')
        fo.write(hdr.rstrip('\n') + '\tlength_range\n')
        for line in f:
            r = dict(zip(h, line.rstrip('\n').split('\t')))
            keep = {k: r[k] for k in ('pair_scope', 'a_sample', 'a_contig', 'a_x_interval', 'b_sample', 'b_contig',
                                      'b_x_interval')}
            wide.append(keep)
            la, lb = int(r['x_len_a']), int(r['x_len_b'])
            if in_target(la) and in_target(lb):
                tgt.append(keep)
                fo.write(line.rstrip('\n') + '\t%d<=len<%d\n' % (TARGET_MIN, TARGET_MAX))
            else:
                below['one_side_below_%d' % TARGET_MIN if (la < TARGET_MIN) != (lb < TARGET_MIN) else
                      'both_sides_below_%d' % TARGET_MIN] += 1
    res = dict(rule='pass both checks and %d <= extracted X length < %d bp on both sides' % (TARGET_MIN, TARGET_MAX),
               table=str(out / name), main=summarize(tgt), wide_lt2500=summarize(wide),
               excluded_from_wide=dict(below))
    (out / ('target_%d_%d_summary.json' % (TARGET_MIN, TARGET_MAX))).write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
