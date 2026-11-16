"""Round-5 cumulative set, ucd_readcheck v0.1 (readcheck_round5/) vs v0.2 (per-interval alleles, this directory).
Same intervals, runs, genome PAFs and thresholds; only the allele alignment differs. Per interval: class, verdict and
count changes; near-threshold intervals of v0.1 (n = P_matched + A_matched in [MIN_SPAN - 1, MIN_SPAN + 2], or the
other allele's share of n within 5 points of MAX_OTHER_FRACTION) and what they became; pair class changes.
Writes interval_changes.tsv and compare.json here. Usage: python compare_v01_v02.py
"""
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import protein_context as pc  # noqa: E402
import scale_entry as se  # noqa: E402,F401
import s3_readcheck as rc3  # noqa: E402

OLD = HERE.parent / 'readcheck_round5'
COUNTS = ('P', 'A', 'P_matched', 'A_matched', 'ambiguous', 'locus_ambiguous', 'off_locus', 'cigar_fail')


def num(x):
    return int(x) if x not in ('', None) else 0


def near_threshold(r):
    """v0.1 interval close to a verdict boundary (expected allele P, as in the round's position check)."""
    p, a = num(r['P_matched']), num(r['A_matched'])
    n = p + a
    if rc3.MIN_SPAN - 1 <= n <= rc3.MIN_SPAN + 2:
        return 'n_near_min_span'
    if n and abs(a / n - rc3.MAX_OTHER_FRACTION) <= 0.05:
        return 'other_share_near_10pct'
    return ''


def main():
    old = {r['test_id']: r for r in pc.read_tsv(OLD / 'positions.tsv.gz')}
    new = {r['test_id']: r for r in pc.read_tsv(HERE / 'positions.tsv.gz')}
    assert set(old) == set(new), (len(set(old) ^ set(new)), 'interval ids differ')
    rows, cls, ver, delta, near = [], Counter(), Counter(), {k: Counter() for k in COUNTS}, Counter()
    for tid in sorted(old):
        o, n = old[tid], new[tid]
        d = {k: num(n[k]) - num(o[k]) for k in COUNTS}
        for k in COUNTS:
            delta[k]['up' if d[k] > 0 else 'down' if d[k] < 0 else 'same'] += 1
        cls[(o['position_class'], n['position_class'])] += 1
        ver[(o['verdict'], n['verdict'])] += 1
        nt = near_threshold(o)
        if nt:
            near[(nt, o['position_class'], n['position_class'])] += 1
        if any(d.values()) or o['position_class'] != n['position_class'] or nt:
            rows.append(dict(test_id=tid, sample=o['sample'], v01_class=o['position_class'], v02_class=n['position_class'],
                             v01_verdict=o['verdict'], v02_verdict=n['verdict'], near_threshold_v01=nt,
                             **{'v01_' + k: o[k] for k in ('P_matched', 'A_matched')},
                             **{'v02_' + k: n[k] for k in ('P_matched', 'A_matched')},
                             **{'d_' + k: d[k] for k in COUNTS}))
    pc.write_tsv(HERE / 'interval_changes.tsv', rows, list(rows[0]))
    op = {r['pair_id']: r['pair_class'] for r in pc.read_tsv(OLD / 'pairs.tsv.gz')}
    np_ = {r['pair_id']: r['pair_class'] for r in pc.read_tsv(HERE / 'pairs.tsv.gz')}
    pcls = Counter((op[k], np_[k]) for k in op)
    sums = {k: (sum(num(old[t][k]) for t in old), sum(num(new[t][k]) for t in new)) for k in COUNTS}
    pm = sorted(num(new[t]['P_matched']) - num(old[t]['P_matched']) for t in old)
    out = dict(intervals=len(old), intervals_with_any_count_change=sum(1 for r in rows if any(r['d_' + k] for k in COUNTS)),
               position_class={'%s -> %s' % k: v for k, v in sorted(cls.items())},
               verdict={'%s -> %s' % k: v for k, v in sorted(ver.items())},
               count_totals_v01_v02=sums, count_change_direction={k: dict(v) for k, v in delta.items()},
               P_matched_change=dict(min=pm[0], max=pm[-1], up=sum(x > 0 for x in pm), down=sum(x < 0 for x in pm)),
               near_threshold_v01={'%s | %s -> %s' % k: v for k, v in sorted(near.items())},
               pair_class={'%s -> %s' % k: v for k, v in sorted(pcls.items())})
    (HERE / 'compare.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
