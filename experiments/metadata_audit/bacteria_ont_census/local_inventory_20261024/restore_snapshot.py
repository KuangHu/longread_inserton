"""Verify and restore the compressed metadata snapshot, without network access."""
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def local(name):
    path = (ROOT / name).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError('Path outside snapshot: ' + name)
    return path

def digest(data):
    return hashlib.sha256(data).hexdigest()

def save(path, data):
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError('Refusing to overwrite different local content: ' + str(path))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)

def main():
    manifest = json.loads((ROOT / 'compressed_manifest.json').read_text(encoding='utf-8'))
    for entry in manifest:
        chunks = []
        for part in entry['parts']:
            compressed = local(part['path']).read_bytes()
            assert len(compressed) == part['bytes'] and digest(compressed) == part['sha256'], part['path']
            chunks.append(gzip.decompress(compressed))
        data = b''.join(chunks)
        assert len(data) == entry['bytes'] and digest(data) == entry['sha256'], entry['output']
        save(local(entry['output']), data)
    cache = local('ena_bacteria_ont_2026-10-24/batch_priority/literature_cache')
    for line in local('current_literature_cache.jsonl').read_text(encoding='utf-8').splitlines():
        item = json.loads(line)
        filename = item['cache_filename']
        if Path(filename).name != filename:
            raise ValueError('Invalid cache filename')
        save(cache / filename, (json.dumps(item['record'],ensure_ascii=False)+'\n').encode('utf-8'))
    print(f'Verified and restored {len(manifest)} metadata files and literature responses.')

if __name__ == '__main__':
    main()
