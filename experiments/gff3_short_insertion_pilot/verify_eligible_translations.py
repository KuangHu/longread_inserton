"""Check eligible CDS translations against the downloaded NCBI protein FASTA."""
import csv
import gzip
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent
details = json.loads((ROOT / 'results/feature_details.json').read_text())
with (ROOT / 'results/candidate_annotations.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
genomes, proteins = {}, {}
for sample in ('H16', 'H850'):
    folder = ROOT / 'data' / sample
    with gzip.open(next(folder.glob('*_genomic.fna.gz')), 'rt') as handle:
        genomes[sample] = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    with gzip.open(next(folder.glob('*_protein.faa.gz')), 'rt') as handle:
        proteins[sample] = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
checks = []
for row in rows:
    if row['annotation_status'] != 'eligible_1_2_CDS':
        continue
    for feature in details[row['seed_id']]:
        if feature['kind'] != 'CDS':
            continue
        sample = row['filled_sample']
        seq = genomes[sample][feature['seqid']].seq
        coding = Seq('').join([seq[a:b] for a, b in sorted(feature['parts'])])
        if feature['strand'] == '-':
            coding = coding.reverse_complement()
        table = int(feature['attrs'].get('transl_table', ['11'])[0])
        translated = str(coding.translate(table=table, cds=True))
        protein = feature['attrs']['protein_id'][0]
        assert translated == str(proteins[sample][protein].seq), (row['seed_id'], protein)
        checks.append({'seed_id': row['seed_id'], 'protein_id': protein,
                       'amino_acids': len(translated), 'translation_matches_FAA': True})
result = {'status': 'pass', 'checked_CDS': len(checks), 'checks': checks,
          'scope': 'Eligible CDS only; no functional or insertion-family validation.'}
(ROOT / 'results/translation_checks.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'status': 'pass', 'checked_CDS': len(checks)}))
