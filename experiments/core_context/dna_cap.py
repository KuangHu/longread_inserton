"""Cap of one family pool for the DNA entry: at most N BioSamples, chosen by a rule that uses assembly metadata only.

Order (fixed before any result is seen; nothing about cores, repeats or reads enters):
  1. group: Complete Genome / Chromosome before Scaffold / Contig;
  2. within a group, round-robin over species (genus + epithet of the organism name; the whole name for "sp."):
     every species gives its 1st assembly, then every species its 2nd, ...; a species' assemblies are ordered by
     level, then accession;
  3. ties by species name, then accession.
The first N rows are kept. NCBI-atypical rows are dropped first. Rows not kept are written with reason "cap".
Reads pool.tsv only; nothing is downloaded here.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

RANK = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}


def species(org):
    t = org.split()
    return org if len(t) < 2 or t[1] in ('sp.', 'cf.', 'aff.') else ' '.join(t[:2])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pool', required=True)
    p.add_argument('--cap', required=True, type=int)
    p.add_argument('--out', required=True, help='capped pool.tsv; <out>.json summary; <out>.dropped.tsv')
    a = p.parse_args()
    rd = csv.DictReader(open(a.pool), delimiter='\t')
    rows = [r for r in rd if not r['atypical']]
    by = defaultdict(list)
    for r in rows:
        r['_g'] = 0 if RANK.get(r['assembly_level'], 9) <= 1 else 1
        by[(r['_g'], species(r['organism']))].append(r)
    for (g, sp), v in by.items():
        v.sort(key=lambda r: (RANK.get(r['assembly_level'], 9), r['download_accession']))
        for i, r in enumerate(v):
            r['_k'] = (g, i, sp, r['download_accession'])
    rows.sort(key=lambda r: r['_k'])
    keep, drop = rows[:a.cap], rows[a.cap:]
    for path, sel in ((a.out, keep), (a.out + '.dropped.tsv', drop)):
        with open(path, 'w', newline='') as f:
            w = csv.DictWriter(f, rd.fieldnames, delimiter='\t', lineterminator='\n', extrasaction='ignore')
            w.writeheader(); w.writerows(sel)
    info = dict(rule='dna_cap_v1', cap=a.cap, pool_not_atypical=len(rows), kept=len(keep), dropped=len(drop),
                kept_levels=dict(Counter(r['assembly_level'] for r in keep)),
                species_in_pool=len({species(r['organism']) for r in rows}),
                species_kept=len({species(r['organism']) for r in keep}))
    Path(a.out + '.json').write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps(info))


if __name__ == '__main__':
    main()
