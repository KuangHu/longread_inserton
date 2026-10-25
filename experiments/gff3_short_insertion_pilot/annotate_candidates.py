"""Small, standalone NCBI GFF3 overlap audit; not an insertion caller.

python annotate_candidates.py --download
Requires Python >=3.10 and Biopython. Inputs/outputs are relative to this file.
Raw NCBI files stay in data/ and must not be committed.
"""
from pathlib import Path
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import platform
from urllib.parse import unquote
from urllib.request import urlopen
import Bio
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent
SOURCES = {
    'H16': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/051/903/345/GCF_051903345.1_ASM5190334v1',
    'H850': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/056/236/025/GCF_056236025.1_ASM5623602v1',
}


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def attributes(value):
    # Split before percent decoding so encoded commas remain literal values.
    return {k: [unquote(v) for v in values.split(',')]
            for item in value.split(';') if '=' in item
            for k, values in [item.split('=', 1)]}


def parse_gff(lines):
    features, regions, directives = {}, {}, []
    for line in lines:
        if line.startswith('##FASTA'):
            break
        if line.startswith('##sequence-region '):
            _, seqid, start, end = line.split()
            regions[seqid] = (int(start), int(end))
        if line.startswith('#!'):
            directives.append(line.strip())
        if not line.strip() or line.startswith('#'):
            continue
        seqid, source, kind, start, end, score, strand, phase, raw = line.rstrip('\n').split('\t')
        attr = attributes(raw)
        if kind not in ('gene', 'pseudogene', 'CDS') and not kind.endswith('RNA'):
            continue
        if 'ID' not in attr:
            raise ValueError(f'Missing feature ID: {seqid}:{start}-{end}')
        key = attr['ID'][0]
        feature = features.setdefault(key, {'id': key, 'seqid': seqid, 'kind': kind,
                                            'strand': strand, 'parts': [], 'attrs': {}})
        if (feature['seqid'], feature['kind'], feature['strand']) != (seqid, kind, strand):
            raise ValueError(f'Inconsistent multipart ID: {key}')
        feature['parts'].append((int(start) - 1, int(end)))
        for k, values in attr.items():
            feature['attrs'].setdefault(k, set()).update(values)
    return features, regions, directives


def overlaps(feature, contig, start, end):
    return feature['seqid'] == contig and any(a < end and b > start for a, b in feature['parts'])


def contained(feature, start, end):
    return all(start <= a < b <= end for a, b in feature['parts'])


def flags(feature, features):
    result, pending, seen = set(), [feature], set()
    while pending:
        f = pending.pop()
        if f['id'] in seen:
            continue
        seen.add(f['id'])
        attrs = f['attrs']
        if f['kind'] == 'pseudogene' or 'true' in attrs.get('pseudo', set()) or 'pseudogene' in attrs:
            result.add('pseudogene')
        if 'true' in attrs.get('partial', set()) or 'start_range' in attrs or 'end_range' in attrs:
            result.add('partial')
        if 'exception' in attrs or 'transl_except' in attrs:
            result.add('translation_exception')
        for parent in attrs.get('Parent', set()):
            if parent not in features:
                result.add('missing_parent')
            else:
                pending.append(features[parent])
    return result


def audit(row, features, proteins):
    start, end = int(row['start0']), int(row['end0'])
    assert end - start == int(row['length'])
    touching = [f for f in features.values() if overlaps(f, row['filled_contig'], start, end)]
    cds = [f for f in touching if f['kind'] == 'CDS']
    genes = [f for f in touching if f['kind'] in ('gene', 'pseudogene')]
    rnas = [f for f in touching if f['kind'].endswith('RNA')]
    whole_cds = [f for f in cds if contained(f, start, end)]
    whole_genes = [f for f in genes if contained(f, start, end)]
    boundary = [f for f in touching if not contained(f, start, end)]
    reasons = set()
    if len(whole_cds) not in (1, 2):
        reasons.add('complete_CDS_count_not_1_or_2')
    if len(whole_genes) != len(whole_cds):
        reasons.add('gene_CDS_count_mismatch')
    if boundary:
        reasons.add('boundary_intersects_feature')
    if rnas:
        reasons.add('RNA_overlap')
    for f in touching:
        reasons.update(flags(f, features))
    protein_ids = sorted({p for f in whole_cds for p in f['attrs'].get('protein_id', [])})
    if any(not f['attrs'].get('protein_id') for f in whole_cds) or any(p not in proteins for p in protein_ids):
        reasons.add('protein_sequence_missing')
    details = []
    for f in sorted(touching, key=lambda f: (min(a for a, b in f['parts']), f['kind'], f['id'])):
        details.append({**f, 'attrs': {k: sorted(v) for k, v in f['attrs'].items()},
                        'contained': contained(f, start, end), 'flags': sorted(flags(f, features))})
    result = {**row, 'complete_gene_count': len(whole_genes), 'complete_CDS_count': len(whole_cds),
              'overlapping_CDS_count': len(cds), 'boundary_feature_count': len(boundary),
              'RNA_feature_count': len(rnas), 'protein_ids': ';'.join(protein_ids),
              'products': ';'.join(sorted({p for f in whole_cds for p in f['attrs'].get('product', [])})),
              'annotation_status': 'review' if reasons else 'eligible_1_2_CDS',
              'review_reasons': ';'.join(sorted(reasons))}
    return result, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    data, output = ROOT / 'data', ROOT / 'results'
    data.mkdir(exist_ok=True)
    output.mkdir(exist_ok=True)
    samples = json.loads((ROOT / 'inputs/samples.json').read_text())
    source_receipts, annotations, protein_maps = [], {}, {}
    for sample in samples:
        name, accession = sample['sample'], sample['assembly']
        base = SOURCES[name]
        prefix = base.rsplit('/', 1)[1]
        names = ['md5checksums.txt'] + [prefix + suffix for suffix in ('_genomic.gff.gz', '_genomic.fna.gz', '_protein.faa.gz')]
        folder = data / name
        folder.mkdir(exist_ok=True)
        files = []
        for filename in names:
            path = folder / filename
            url = base + '/' + filename
            if args.download and not path.exists():
                temporary = path.with_suffix(path.suffix + '.partial')
                with urlopen(url, timeout=120) as response, temporary.open('wb') as handle:
                    while chunk := response.read(1024 * 1024):
                        handle.write(chunk)
                temporary.replace(path)
            if not path.exists():
                raise FileNotFoundError(f'{path}: rerun with --download')
            files.append({'filename': filename, 'url': url, 'bytes': path.stat().st_size, 'sha256': digest(path)})
        md5s = {name.removeprefix('./'): md5 for md5, name in
                (line.split() for line in (folder / names[0]).read_text().splitlines() if line.strip())}
        for receipt in files[1:]:
            assert digest(folder / receipt['filename'], 'md5') == md5s[receipt['filename']]
            receipt['archive_md5'] = md5s[receipt['filename']]
        with gzip.open(folder / names[2], 'rt') as handle:
            records = {r.id: str(r.seq).upper() for r in SeqIO.parse(handle, 'fasta')}
        assert set(records) == {c['accession'] for c in sample['contigs']}
        for c in sample['contigs']:
            seq = records[c['accession']]
            assert len(seq) == c['length'] and hashlib.sha256(seq.encode()).hexdigest() == c['sha256']
        with gzip.open(folder / names[1], 'rt') as handle:
            features, regions, directives = parse_gff(handle)
        assert any(accession in d for d in directives if d.startswith('#!genome-build-accession'))
        assert regions == {c['accession']: (1, c['length']) for c in sample['contigs']}
        assert all(f['seqid'] in records for f in features.values())
        for f in features.values():
            assert all(0 <= a < b <= len(records[f['seqid']]) for a, b in f['parts']), 'Origin-spanning/out-of-range feature requires separate handling'
        with gzip.open(folder / names[3], 'rt') as handle:
            proteins = {r.id: str(r.seq) for r in SeqIO.parse(handle, 'fasta')}
        annotations[name], protein_maps[name] = features, proteins
        source_receipts.append({'sample': name, 'assembly': accession, 'files': files,
                                'directives': directives, 'sequence_sha256_matches_frozen': True,
                                'feature_counts': dict(Counter(f['kind'] for f in features.values())),
                                'protein_fasta_records': len(proteins)})
    with (ROOT / 'inputs/assembly_candidates.tsv').open() as handle:
        candidates = list(csv.DictReader(handle, delimiter='\t'))
    selected = [r for r in candidates if 0 < int(r['length']) < 2500]
    rows, details = [], {}
    for row in selected:
        annotated, features = audit(row, annotations[row['filled_sample']], protein_maps[row['filled_sample']])
        rows.append(annotated)
        details[row['seed_id']] = features
    with (output / 'candidate_annotations.tsv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    summary = {'input_candidates': len(candidates), 'core_length_lt_2500': len(rows),
               'status_counts': dict(Counter(r['annotation_status'] for r in rows)),
               'complete_CDS_count_distribution': dict(Counter(r['complete_CDS_count'] for r in rows)),
               'eligible_CDS_count_distribution': dict(Counter(r['complete_CDS_count'] for r in rows if r['annotation_status'] == 'eligible_1_2_CDS')),
               'review_reason_counts_nonexclusive': dict(Counter(reason for r in rows for reason in r['review_reasons'].split(';') if reason)),
               'limits': ['Input discovery used 400-20000 bp; candidates below 400 bp were not searched.',
                          'Annotation eligibility is not insertion validation, family count or independent event count.',
                          'Current annotation snapshot; annotation may change without assembly version change.',
                          'No new read validation or protein clustering in this audit.']}
    for filename, value in [('summary.json', summary), ('feature_details.json', details),
                            ('provenance.json', {'retrieved_or_checked_utc': datetime.now(timezone.utc).isoformat(),
                                                'input_sha256': {p.name: digest(p) for p in (ROOT / 'inputs').iterdir() if p.is_file()},
                                                'script_sha256': digest(__file__),
                                                'python_version': platform.python_version(),
                                                'biopython_version': Bio.__version__, 'sources': source_receipts})]:
        (output / filename).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
