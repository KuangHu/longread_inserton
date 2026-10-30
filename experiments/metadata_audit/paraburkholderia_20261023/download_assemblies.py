"""Download matched GCA genomic FASTA packages for this metadata/ANI audit only."""
import hashlib, io, json, urllib.request, zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'genomes'
OUT.mkdir(exist_ok=True)

def fetch(g):
    acc = g['accession']
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA'
    result = {'accession': acc, 'url': url}
    for attempt in range(1, 3):
        try:
            request = urllib.request.Request(url, headers={'Accept': 'application/zip'})
            data = urllib.request.urlopen(request, timeout=60).read()
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                bad = z.testzip()
                if bad:
                    raise ValueError('ZIP CRC mismatch: ' + bad)
                files = [n for n in z.namelist() if n.endswith('_genomic.fna')]
                if len(files) != 1 or acc not in files[0]:
                    raise ValueError('Unexpected FASTA members: ' + str(files))
                seq = z.read(files[0])
            target = OUT / (acc + '.fna')
            target.write_bytes(seq)
            lines = seq.splitlines()
            bases = sum(len(x.strip()) for x in lines if not x.startswith(b'>'))
            expected = int(g['assembly_stats']['total_sequence_length'])
            if bases != expected:
                raise ValueError(f'Length mismatch: {bases} != {expected}')
            result.update(status='downloaded_verified_length', attempt=attempt,
                          fasta_bytes=len(seq), bases=bases, records=sum(x.startswith(b'>') for x in lines),
                          fasta_sha256=hashlib.sha256(seq).hexdigest(), zip_sha256=hashlib.sha256(data).hexdigest())
            return result
        except Exception as exc:
            result.update(status='download_failed', attempt=attempt, error=str(exc))
    return result

if __name__ == '__main__':
    gs = {g['accession']:g for g in json.loads((ROOT/'matched_assemblies.json').read_text())}
    results = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(fetch,g) for g in gs.values()]):
            r = future.result()
            results.append(r)
            print(json.dumps(r), flush=True)
            (ROOT/'assembly_download_audit.json').write_text(json.dumps(sorted(results,key=lambda x:x['accession']),indent=2)+'\n')
