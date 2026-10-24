"""Rename archive FASTQ records to `<molecule>/<index>` and diagnose the read unit.

`rename` streams FASTQ to stdout so minimap2 query names carry the molecule:
PacBio subreads keep movie/ZMW, SRA spots are used when the original name was
stripped (one spot is one ZMW for PacBio RS loads), and each Nanopore read is
its own molecule. `diagnose` scans records for SMRTbell adapter sequence, which
distinguishes whole polymerase reads from adapter-split subreads.
"""

import argparse
import gzip
import json
import re
import sys

SUBREAD = re.compile(r"^(m[^/\s]+/\d+)/(?:\d+_\d+|ccs)$")
# SMRTbell hairpin adapter core; checked on both strands.
ADAPTER = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"
COMP = str.maketrans("ACGT", "TGCA")


def molecule_of(header):
    fields = header[1:].split()
    spot = fields[0]
    if len(fields) > 1:
        match = SUBREAD.match(fields[1])
        if match:
            return match.group(1)
    return spot


def records(path):
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().rstrip("\n")
            handle.readline()
            qual = handle.readline().rstrip("\n")
            yield header.rstrip("\n"), seq, qual


def adapter_regions(seq, kmers, k=15, gap=60):
    """Clusters of >=2 distinct adapter 15-mers, as (start, end) on the read."""
    hits = []
    for kmer in kmers:
        i = seq.find(kmer)
        while i >= 0:
            hits.append((i, kmer))
            i = seq.find(kmer, i + 1)
    hits.sort()
    regions, cluster = [], []
    for pos, kmer in hits + [(None, None)]:
        if cluster and (pos is None or pos - cluster[-1][0] > gap):
            if len({km for _, km in cluster}) >= 2:
                regions.append((cluster[0][0], cluster[-1][0] + k))
            cluster = []
        if pos is not None:
            cluster.append((pos, kmer))
    return regions


def rename(path, split_adapters, min_piece, stats_path):
    kmers = adapter_kmers()
    counts = {}
    stats = dict(records=0, records_split=0, adapter_regions=0, pieces_written=0, pieces_dropped_short=0)
    out = sys.stdout
    for header, seq, qual in records(path):
        stats["records"] += 1
        mol = molecule_of(header).replace("/", ":")
        pieces = [(0, len(seq))]
        if split_adapters:
            regions = adapter_regions(seq, kmers)
            if regions:
                stats["records_split"] += 1
                stats["adapter_regions"] += len(regions)
                bounds = [0] + [x for r in regions for x in r] + [len(seq)]
                pieces = list(zip(bounds[::2], bounds[1::2]))
        for start, end in pieces:
            if end - start < min_piece:
                stats["pieces_dropped_short"] += 1
                continue
            counts[mol] = counts.get(mol, 0) + 1
            stats["pieces_written"] += 1
            out.write("@%s/%d\n%s\n+\n%s\n" % (mol, counts[mol], seq[start:end], qual[start:end]))
    stats["molecules"] = len(counts)
    if stats_path:
        with open(stats_path, "w") as handle:
            json.dump(stats, handle, indent=2)


def adapter_kmers(k=15):
    rc = ADAPTER.translate(COMP)[::-1]
    return {seq[i : i + k] for seq in (ADAPTER, rc) for i in range(len(seq) - k + 1)}


def diagnose(path, limit):
    # At CLR accuracy a full adapter is rarely error-free, but several of its
    # 15-mers usually are; two distinct hits is far above random background.
    kmers = adapter_kmers()
    n = with_adapter = 0
    lengths = []
    molecules = set()
    for header, seq, _ in records(path):
        n += 1
        molecules.add(molecule_of(header))
        lengths.append(len(seq))
        with_adapter += sum(1 for k in kmers if k in seq) >= 2
        if n >= limit:
            break
    lengths.sort()
    fraction = with_adapter / n if n else 0
    return dict(
        records_scanned=n,
        molecules_in_scan=len(molecules),
        records_with_adapter=with_adapter,
        fraction_with_adapter=round(fraction, 4),
        rule=">=2 distinct SMRTbell adapter 15-mers (either strand) in the record",
        median_length=lengths[len(lengths) // 2] if lengths else None,
        interpretation=(
            "adapters_inside_records: likely unsplit polymerase reads"
            if fraction > 0.05
            else "few_internal_adapters: adapter-split subreads or adapter-free reads"
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["rename", "diagnose"])
    parser.add_argument("fastq")
    parser.add_argument("--limit", type=int, default=20000)
    parser.add_argument("--split-adapters", action="store_true", help="Cut PacBio records at internal SMRTbell adapters")
    parser.add_argument("--min-piece", type=int, default=500)
    parser.add_argument("--stats", help="Write rename statistics JSON here")
    args = parser.parse_args()
    if args.mode == "rename":
        rename(args.fastq, args.split_adapters, args.min_piece, args.stats)
    else:
        print(json.dumps(diagnose(args.fastq, args.limit), indent=2))


if __name__ == "__main__":
    main()
