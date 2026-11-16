"""DNA-entry manifest of one downloaded family pool (fna_pool/build_pool.py + download_fna.py), without the S2 ANI step.

One row per pool BioSample whose download accession is verified (verify.json status verified, every check true) and
not NCBI-atypical. Columns are those of scale_joint_v0.2/manifest_b1-8r.tsv; batch is 1 for every row;
species_cluster is the organism name and clone_group the BioSample (no ANI clustering was run: the DNA entry uses
species_cluster only to order assemblies into parts). GFF3 / FAA columns are filled when the files exist.
Reads local files only; nothing is downloaded here.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

COLS = ['batch', 'sample', 'assembly', 'assembly_level', 'genus', 'organism', 'species_cluster', 'clone_group', 'fasta',
        'fna_sha256', 'gff3', 'gff3_sha256', 'faa', 'faa_sha256', 'total_bp', 'records', 'family', 'level_batch']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pool', required=True, help='pool.tsv of build_pool.py build')
    p.add_argument('--data', required=True, type=Path, help='data root of download_fna.py')
    p.add_argument('--family', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    rows, skipped = [], Counter()
    for r in csv.DictReader(open(a.pool), delimiter='\t'):
        if r['atypical']:
            skipped['atypical'] += 1; continue
        d = a.data / 'genomes' / r['download_accession']
        if not (d / 'verify.json').exists():
            skipped['not_downloaded'] += 1; continue
        v = json.loads((d / 'verify.json').read_text())
        if v['status'] != 'verified' or not all(v['checks'].values()) or v['biosample'] != r['biosample']:
            skipped['verify_failed'] += 1; continue
        sha = {f['file']: f['sha256'] for f in v['files']}
        fna = [f for f in sha if f.endswith('_genomic.fna')]
        if len(fna) != 1 or not (d / fna[0]).is_file():
            skipped['no_fna'] += 1; continue
        opt = lambda f: (str(d / f), sha[f]) if f in sha and (d / f).is_file() else ('', '')
        gff, faa = opt('genomic.gff'), opt('protein.faa')
        rows.append(dict(batch=1, sample=r['biosample'], assembly=r['download_accession'], assembly_level=r['assembly_level'],
                         genus=r['organism'].split()[0], organism=r['organism'], species_cluster=r['organism'],
                         clone_group=r['biosample'], fasta=str(d / fna[0]), fna_sha256=sha[fna[0]], gff3=gff[0],
                         gff3_sha256=gff[1], faa=faa[0], faa_sha256=faa[1], total_bp=v['total_bp'], records=v['records'],
                         family=a.family, level_batch=''))
    assert len({m['sample'] for m in rows}) == len(rows), 'a BioSample twice in the pool'
    rows.sort(key=lambda m: m['assembly'])
    with open(a.out + '.tmp', 'w', newline='') as f:
        w = csv.DictWriter(f, COLS, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    Path(a.out + '.tmp').replace(a.out)
    info = dict(family=a.family, assemblies=len(rows), skipped=dict(skipped),
                levels=dict(Counter(m['assembly_level'] for m in rows)))
    Path(a.out + '.json').write_text(json.dumps(info, indent=1) + '\n')
    print(json.dumps(info))


if __name__ == '__main__':
    main()
