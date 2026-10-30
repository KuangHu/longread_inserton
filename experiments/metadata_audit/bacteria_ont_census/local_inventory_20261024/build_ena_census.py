"""Snapshot public ENA metadata, without downloading sequence files or selecting samples."""
import argparse
import collections
import hashlib
import json
import re
import sqlite3
import subprocess
import shutil
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, urlsplit

ROOT = Path(__file__).resolve().parent / 'ena_bacteria_ont_2026-10-24'
QUERY = 'tax_tree(2) AND instrument_platform="OXFORD_NANOPORE"'
FIELDS = ['run_accession', 'experiment_accession', 'sample_accession',
          'secondary_sample_accession', 'study_accession', 'secondary_study_accession',
          'tax_id', 'scientific_name', 'instrument_platform', 'instrument_model',
          'library_strategy', 'library_source', 'library_selection', 'library_layout',
          'read_count', 'base_count', 'first_public', 'last_updated',
          'fastq_ftp', 'fastq_bytes', 'fastq_md5', 'submitted_ftp',
          'submitted_bytes', 'submitted_md5', 'submitted_format']


def fetch(endpoint, params, target, seconds=900):
    url = 'https://www.ebi.ac.uk/ena/portal/api/' + endpoint + '?' + urlencode(params)
    tmp = target.with_suffix(target.suffix + '.partial')
    print('Fetching metadata:', target.name, flush=True)
    subprocess.run([shutil.which('curl.exe') or shutil.which('curl') or 'curl', '--silent', '--show-error', '--fail', '--location',
                    '--retry', '2', '--retry-delay', '3', '--connect-timeout', '20',
                    '--max-time', str(seconds), '--max-filesize', '1073741824',
                    '--output', str(tmp), url], check=True)
    json.loads(tmp.read_text(encoding='utf-8'))
    tmp.replace(target)
    return url


def download():
    ROOT.mkdir(exist_ok=True)
    start = datetime.now(timezone.utc).isoformat()
    params = dict(result='read_run', query=QUERY, format='json')
    before = fetch('count', params, ROOT / 'count_before.json', 60)
    url = fetch('search', dict(**params, fields=','.join(FIELDS), limit=0), ROOT / 'ena_raw.json')
    after = fetch('count', params, ROOT / 'count_after.json', 60)
    provenance = dict(started_utc=start, completed_utc=datetime.now(timezone.utc).isoformat(),
                      query=QUERY, source_url=url, count_urls=[before, after], fields=FIELDS,
                      sequence_files_downloaded=0,
                      scope='Public ENA read_run records under Bacteria taxonomy, all ONT library types; excludes records classified outside tax_tree(2).')
    (ROOT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')


def split(value):
    return value.split(';') if value else []


def file_entries(row, kind):
    urls, sizes, digests = [split(row.get(kind + suffix, '')) for suffix in ['_ftp', '_bytes', '_md5']]
    aligned = len(urls) == len(sizes) == len(digests)
    out = []
    for i, url in enumerate(urls):
        size = sizes[i] if i < len(sizes) else ''
        digest = digests[i] if i < len(digests) else ''
        normalized = 'https://' + url if '://' not in url else url
        parsed = urlsplit(normalized)
        out.append(dict(run_accession=row['run_accession'], file_kind=kind, file_index=i + 1,
                        source_url=url, download_url=normalized,
                        expected_bytes=int(size) if size.isdigit() else None,
                        expected_md5=digest,
                        metadata_complete=aligned and size.isdigit() and int(size) > 0 and bool(re.fullmatch('[0-9a-fA-F]{32}', digest)) and bool(parsed.hostname),
                        remote_file_status='not_probed', local_file_status='not_checked'))
    return out, aligned


def build():
    raw = ROOT / 'ena_raw.json'
    rows = json.loads(raw.read_text(encoding='utf-8'))
    expected = [int(json.loads((ROOT / ('count_' + x + '.json')).read_text())['count']) for x in ['before', 'after']]
    ids = [r['run_accession'] for r in rows]
    assert len(ids) == len(set(ids)), 'Duplicate runs in API export'
    assert all(r['instrument_platform'] == 'OXFORD_NANOPORE' for r in rows)
    assert all(re.fullmatch('[SED]RR[0-9]+', r) for r in ids)
    old_path = ROOT.parent / 'run_catalog.json'
    old_data = json.loads(old_path.read_text(encoding='utf-8')) if old_path.exists() else {}
    old_rows = old_data if isinstance(old_data, list) else old_data.get('runs', [])
    old_ids = {r.get('run', r.get('run_accession', '')) for r in old_rows}
    db_path = ROOT / 'run_inventory.sqlite'
    conn = sqlite3.connect(db_path)
    conn.executescript('DROP TABLE IF EXISTS runs; DROP TABLE IF EXISTS files;')
    extra = ['availability', 'fastq_file_count', 'fastq_total_bytes', 'submitted_file_count',
             'fastq_metadata_complete', 'submitted_metadata_complete', 'assembly_link_status',
             'in_previous_registry', 'metadata_json']
    conn.execute('CREATE TABLE runs (' + ','.join('"' + f + '" TEXT' for f in FIELDS + extra) + ', PRIMARY KEY(run_accession))')
    conn.execute('CREATE TABLE files (run_accession TEXT, file_kind TEXT, file_index INTEGER, source_url TEXT, download_url TEXT, expected_bytes INTEGER, expected_md5 TEXT, metadata_complete INTEGER, remote_file_status TEXT, local_file_status TEXT, PRIMARY KEY(run_accession,file_kind,file_index))')
    statuses = collections.Counter()
    strategies = collections.Counter()
    source_counts = collections.Counter()
    models = collections.Counter()
    totals = collections.Counter()
    samples, studies = set(), set()
    issues = []
    with (ROOT / 'all_runs.jsonl').open('w', encoding='utf-8') as all_out, \
         (ROOT / 'fastq_available_runs.jsonl').open('w', encoding='utf-8') as fastq_out, \
         (ROOT / 'file_inventory.jsonl').open('w', encoding='utf-8') as file_out:
        for r in sorted(rows, key=lambda r: r['run_accession']):
            fastq, fa = file_entries(r, 'fastq')
            submitted, sa = file_entries(r, 'submitted')
            status = 'fastq_listed' if fastq else 'submitted_only' if submitted else 'no_file_link_listed'
            fq_ok = bool(fastq) and all(f['metadata_complete'] for f in fastq)
            sub_ok = bool(submitted) and all(f['metadata_complete'] for f in submitted)
            fq_size = sum(f['expected_bytes'] for f in fastq) if fastq and all(f['expected_bytes'] is not None for f in fastq) else None
            derived = dict(availability=status, fastq_file_count=len(fastq), fastq_total_bytes=fq_size,
                           submitted_file_count=len(submitted), fastq_metadata_complete=fq_ok,
                           submitted_metadata_complete=sub_ok, assembly_link_status='not_assessed',
                           in_previous_registry=r['run_accession'] in old_ids)
            record = dict(**r, **derived)
            text = json.dumps(record, ensure_ascii=False)
            all_out.write(text + '\n')
            if fastq:
                fastq_out.write(text + '\n')
            vals = [r.get(f, '') for f in FIELDS] + [derived[k] for k in extra[:-1]] + [text]
            conn.execute('INSERT INTO runs VALUES (' + ','.join('?' for _ in vals) + ')', vals)
            for f in fastq + submitted:
                file_out.write(json.dumps(f, ensure_ascii=False) + '\n')
                conn.execute('INSERT INTO files VALUES (?,?,?,?,?,?,?,?,?,?)', list(f.values()))
            statuses[status] += 1
            strategies[r.get('library_strategy', '')] += 1
            source_counts[r.get('library_source', '')] += 1
            models[r.get('instrument_model', '')] += 1
            if r.get('sample_accession'): samples.add(r['sample_accession'])
            if r.get('study_accession'): studies.add(r['study_accession'])
            totals['fastq_files'] += len(fastq)
            totals['submitted_files'] += len(submitted)
            totals['fastq_runs_complete_metadata'] += fq_ok
            totals['known_fastq_bytes_sum_per_run'] += fq_size or 0
            if not fa or not sa or (fastq and not fq_ok) or (submitted and not sub_ok):
                issues.append(dict(run_accession=r['run_accession'], fastq_arrays_aligned=fa, submitted_arrays_aligned=sa,
                                   fastq_metadata_complete=fq_ok, submitted_metadata_complete=sub_ok))
    conn.executescript('CREATE INDEX runs_sample ON runs(sample_accession); CREATE INDEX runs_availability ON runs(availability); CREATE INDEX runs_strategy ON runs(library_strategy); CREATE INDEX files_url ON files(download_url);')
    conn.commit()
    dedup = conn.execute("SELECT COUNT(*), SUM(expected_bytes) FROM (SELECT download_url, MAX(expected_bytes) expected_bytes FROM files WHERE file_kind='fastq' GROUP BY download_url)").fetchone()
    id_lists = {
        'all_run_accessions.txt': '1=1',
        'runs_with_file_links.txt': "availability != 'no_file_link_listed'",
        'fastq_run_accessions.txt': "availability = 'fastq_listed'",
        'submitted_only_run_accessions.txt': "availability = 'submitted_only'",
        'no_file_link_run_accessions.txt': "availability = 'no_file_link_listed'",
    }
    for name, condition in id_lists.items():
        records = conn.execute('SELECT run_accession FROM runs WHERE ' + condition + ' ORDER BY run_accession')
        with (ROOT / name).open('w', encoding='utf-8', newline='\n') as out:
            for (accession,) in records:
                out.write(accession + '\n')
    wgs_fastq = conn.execute("SELECT COUNT(*) FROM runs WHERE library_strategy='WGS' AND availability='fastq_listed'").fetchone()[0]
    integrity = conn.execute('PRAGMA integrity_check').fetchone()[0]
    conn.close()
    counts_match = len(rows) == expected[0] == expected[1]
    summary = dict(run_count=len(rows), count_before=expected[0], count_after=expected[1], count_reconciliation_passed=counts_match,
                   unique_runs=len(set(ids)), unique_sample_accessions=len(samples), unique_studies=len(studies),
                   availability=dict(statuses), library_strategy=dict(strategies), library_source=dict(source_counts),
                   instrument_model=dict(models), totals=dict(totals), distinct_fastq_urls=dedup[0],
                   unique_url_known_fastq_bytes=dedup[1], metadata_issue_runs=len(issues),
                   sqlite_integrity=integrity, raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),
                   previous_registry_overlap=len(old_ids.intersection(ids)), assembly_link_status='not_assessed',
                   wgs_fastq_runs=wgs_fastq,
                   remote_file_validation='URLs listed in ENA metadata; not individually probed',
                   scope=QUERY)
    (ROOT / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    (ROOT / 'metadata_issues.json').write_text(json.dumps(issues, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
    assert counts_match, 'API counts changed or export incomplete; see summary'
    assert integrity == 'ok'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fetch', action='store_true')
    args = parser.parse_args()
    if args.fetch:
        download()
    build()
