"""Compare the check results of the regression variants (counts per test_id) and with rounds 4 and 5."""
import csv
import sys
from pathlib import Path

D = Path(sys.argv[1])
R = Path(sys.argv[2])  # readcheck dir holding round4/, round5/
S = 'SAMN31258398'
COLS = ['verdict', 'P', 'A', 'P_matched', 'A_matched', 'ambiguous', 'locus_ambiguous', 'off_locus', 'cigar_fail']


def load(path):
    return {r['test_id']: tuple(r[c] for c in COLS) for r in csv.DictReader(open(str(path)), delimiter='\t')}


res = {v: load(D / v / 'check' / S / 'positions_result.tsv') for v in (D / 'variants.txt').read_text().split()}
res['round4'] = load(R / 'round4' / 'check' / S / 'positions_result.tsv')
res['round5'] = load(R / 'round5' / 'check' / S / 'positions_result.tsv')
base = res['full_a']
out = []
for v, r in res.items():
    diff = sorted(t for t in r if r[t] != base.get(t))
    out.append(dict(variant=v, intervals=len(r), differ_from_full_a=len(diff), test_ids=';'.join(diff)))
    print(v, len(r), 'differ from full_a:', len(diff), diff[:6])
for t in ('P_1ddca118bf21', 'P_1bd7640586ad', 'P_28b0e2c768be'):
    print(t)
    for v, r in res.items():
        if t in r:
            print('   %-22s %s' % (v, dict(zip(COLS, r[t]))))
with open(str(D / 'compare.tsv'), 'w') as h:
    w = csv.DictWriter(h, list(out[0]), delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(out)
