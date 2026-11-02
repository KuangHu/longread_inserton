"""Build the round directories of the SAMN31258398 order / allele-set regression (one positions.tsv each)."""
import csv
import random
import sys
from pathlib import Path

src, out = Path(sys.argv[1]), Path(sys.argv[2])  # round-5 positions.tsv, output directory
S = 'SAMN31258398'
rows = [r for r in csv.DictReader(open(str(src)), delimiter='\t') if r['sample'] == S]
cols = list(rows[0])
cluster = [r for r in rows if r['contig'] == 'NZ_CP109906.1' and 5195190 <= int(r['lo']) <= 5195215]
variants = dict(full_a=rows, full_b=rows, full_shuffle1=rows[:], full_shuffle2=rows[:], cluster=cluster)
random.Random(1).shuffle(variants['full_shuffle1'])
random.Random(2).shuffle(variants['full_shuffle2'])
for t in ('P_1ddca118bf21', 'P_1bd7640586ad', 'P_28b0e2c768be'):
    variants['alone_' + t] = [r for r in rows if r['test_id'] == t]
for name, rs in variants.items():
    d = out / name
    d.mkdir(parents=True, exist_ok=True)
    with open(str(d / 'positions.tsv'), 'w') as h:
        w = csv.DictWriter(h, cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rs)
    print(name, len(rs))
(out / 'variants.txt').write_text('\n'.join(variants) + '\n')
