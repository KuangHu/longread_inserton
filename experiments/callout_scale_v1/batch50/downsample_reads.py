"""Deterministic per-molecule downsample of a verified ONT FASTQ to a target depth, for mapping only.

Each FASTQ record is one molecule (verify.json molecule_grouping must say so). A record's key is
SHA256("<seed>:<read name>"); records are taken in key order until their bases reach
target_depth x assembly bp, so the same seed, input and target always give the same reads,
independent of file order. The full FASTQ stays untouched; discovery never uses reads.
Writes <out>.fastq.gz and <out>.json (seed, threshold key, input/output counts and SHA256).
"""

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def records(path):
    with gzip.open(path, "rb") as h:
        while True:
            head = h.readline()
            if not head:
                return
            seq, plus, qual = h.readline(), h.readline(), h.readline()
            if not head.startswith(b"@") or not plus.startswith(b"+") or not qual:
                raise SystemExit("Malformed FASTQ near " + head[:80].decode(errors="replace"))
            yield head, seq, plus, qual


def key(seed, head):
    return hashlib.sha256(b"%d:" % seed + head[1:].split()[0]).digest()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verify", type=Path, required=True, help="fetch_long_reads.py verify.json of the full FASTQ")
    p.add_argument("--assembly-bp", type=int, required=True)
    p.add_argument("--depth", type=float, default=80.0)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out", type=Path, required=True, help="output prefix")
    a = p.parse_args()
    v = json.loads(a.verify.read_text())
    if v["status"] != "verified" or v["summary"]["molecule_grouping"] != "one_record_per_molecule":
        raise SystemExit("Input must be a verified FASTQ with one record per molecule")
    src = Path(v["local_file"])
    target = int(a.depth * a.assembly_bp)

    keyed, names = [], set()
    for head, seq, _, _ in records(src):
        name = head[1:].split()[0]
        if name in names:
            raise SystemExit("Duplicate read name " + name.decode())
        names.add(name)
        keyed.append((key(a.seed, head), len(seq.rstrip(b"\r\n"))))
    names = None
    keyed.sort()
    taken, bases, threshold = 0, 0, None
    for k, n in keyed:
        if bases >= target:
            break
        taken, bases, threshold = taken + 1, bases + n, k
    in_reads, in_bases = len(keyed), sum(n for _, n in keyed)
    keyed = None

    fq = Path(str(a.out) + ".fastq.gz")
    out_reads = out_bases = 0
    # Fixed gzip header (no name, mtime 0) so the same selection always has the same SHA256.
    with open(fq, "wb") as raw, gzip.GzipFile(filename="", mode="wb", compresslevel=1, fileobj=raw, mtime=0) as w:
        for rec in records(src):
            if threshold is not None and key(a.seed, rec[0]) <= threshold:
                w.writelines(rec)
                out_reads += 1
                out_bases += len(rec[1].rstrip(b"\r\n"))
    if (out_reads, out_bases) != (taken, bases):
        raise SystemExit("Second pass selected a different set")
    report = dict(method="sha256(seed:read_name) order, take until target bases", seed=a.seed,
                  target_depth=a.depth, assembly_bp=a.assembly_bp, target_bases=target,
                  threshold_key=threshold.hex() if threshold else None,
                  input_file=str(src), input_sha256=v["sha256"], input_reads=in_reads, input_bases=in_bases,
                  input_depth=round(in_bases / a.assembly_bp, 1), output_file=str(fq), output_reads=out_reads,
                  output_bases=out_bases, output_depth=round(out_bases / a.assembly_bp, 1),
                  output_sha256=sha256(fq), downsampled=out_reads < in_reads)
    Path(str(a.out) + ".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
