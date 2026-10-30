"""Pick small comparable assembly subsets for the graph entry from measured skani distances.

Inputs are the verified fna_pool stage-1 Comamonadaceae catalog (S2 catalog/groups), its skani
triangle edges and the per-genome verify.json; nothing is downloaded. Only Complete Genome,
non-atypical assemblies with verified FNA + GFF3 + FAA are eligible. For each seed species cluster
the anchor is the member with the most near (ANI >= 99) eligible partners from other BioProjects,
ties broken by accession; the other members are taken per ANI-to-anchor bin,
each from a BioProject and clone group not yet in the subset, in accession order. A bin that cannot
be filled is reported as a shortfall, never back-filled from another bin.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

BINS = [('near', 99.0, 100.01, 2), ('intermediate', 95.0, 99.0, 2), ('medium', 80.0, 95.0, 2)]


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', type=Path, required=True)
    p.add_argument('--edges', type=Path, required=True)
    p.add_argument('--genomes', type=Path, required=True, help='data/genomes directory of the pool')
    p.add_argument('--seed', action='append', required=True, help='name=species_cluster')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    rows = list(csv.DictReader(open(a.catalog), delimiter='\t'))
    eligible, reasons = {}, defaultdict(int)
    for r in rows:
        d = a.genomes / r['assembly']
        if r['assembly_level'] != 'Complete Genome' or r['atypical']:
            reasons['not_complete_or_atypical'] += 1; continue
        if not (d / 'genomic.gff').exists() or not (d / 'protein.faa').exists():
            reasons['no_gff3_faa'] += 1; continue
        v = json.loads((d / 'verify.json').read_text())
        if v['status'] != 'verified' or not all(v['checks'].values()):
            reasons['verify_failed'] += 1; continue
        eligible[r['assembly']] = r
    ani = defaultdict(dict)
    with open(a.edges) as f:
        for e in csv.DictReader(f, delimiter='\t'):
            x, y = (Path(e[k]).parent.name for k in ('Ref_file', 'Query_file'))
            if x in eligible and y in eligible:
                v = (float(e['ANI']), float(e['Align_fraction_ref']), float(e['Align_fraction_query']))
                ani[x][y] = v; ani[y][x] = (v[0], v[2], v[1])
    manifest, notes = [], []
    for seed in a.seed:
        name, cluster = seed.split('=')
        members = sorted(x for x, r in eligible.items() if r['species_cluster'] == cluster)
        near = lambda x: len({eligible[y]['bioproject'] for y, v in ani[x].items()
                              if v[0] >= 99.0 and eligible[y]['bioproject'] != eligible[x]['bioproject']})
        anchor = min(members, key=lambda x: (-near(x), x))
        genus = eligible[anchor]['genus']
        chosen = [(anchor, 'anchor', 100.0, 100.0, 100.0)]
        projects, clones = {eligible[anchor]['bioproject']}, {eligible[anchor]['clone_group']}
        for label, lo, hi, n in BINS:
            pool = sorted(y for y, v in ani[anchor].items() if lo <= v[0] < hi and eligible[y]['genus'] == genus)
            got = 0
            for y in pool:
                r = eligible[y]
                if got == n or r['bioproject'] in projects or r['clone_group'] in clones:
                    continue
                chosen.append((y, label) + ani[anchor][y]); projects.add(r['bioproject']); clones.add(r['clone_group']); got += 1
            notes.append(dict(subset=name, bin=label, ani_range=[lo, min(hi, 100.0)], wanted=n, chosen=got,
                              candidates_in_range=len(pool), shortfall=n - got))
        for acc, label, v, afa, afq in chosen:
            r, d = eligible[acc], a.genomes / acc
            fna = next(d.glob('*_genomic.fna'))
            manifest.append(dict(subset=name, sample=r['sample_id'], assembly=acc, fasta=str(fna), fna_sha256=digest(fna),
                gff3=str(d / 'genomic.gff'), gff3_sha256=digest(d / 'genomic.gff'), faa=str(d / 'protein.faa'),
                faa_sha256=digest(d / 'protein.faa'), organism=r['organism'], species_cluster=r['species_cluster'],
                clone_group=r['clone_group'], bioproject=r['bioproject'], records=r['records'], total_bp=r['total_bp'],
                anchor=anchor, distance_bin=label, ani_to_anchor=v, af_anchor=afa, af_member=afq))
        if fna_mismatch := [m['assembly'] for m in manifest if m['subset'] == name
                            and m['fna_sha256'] != eligible[m['assembly']]['fna_sha256']]:
            raise ValueError('FNA hash differs from catalog: ' + ','.join(fna_mismatch))
    cols = list(manifest[0])
    for name in {m['subset'] for m in manifest}:
        with open(a.out / f'{name}.samples.tsv', 'w', newline='') as f:
            w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n'); w.writeheader()
            w.writerows(m for m in manifest if m['subset'] == name)
    (a.out / 'selection.json').write_text(json.dumps(dict(
        rule=__doc__.strip(), bins=BINS, catalog=str(a.catalog), catalog_sha256=digest(a.catalog),
        edges=str(a.edges), edges_sha256=digest(a.edges), seeds=a.seed, eligible=len(eligible),
        not_eligible=dict(reasons), bins_filled=notes, downloads=0), indent=2) + '\n')


if __name__ == '__main__':
    main()
