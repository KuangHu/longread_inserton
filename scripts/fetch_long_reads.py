"""Download and verify the registered long-read runs (docs/LONG_READ_AVAILABILITY.md).

`probe` reads only the first bytes of each archive FASTQ to classify read names,
read unit and quality encoding. `fetch` downloads, checks archive size/MD5,
records SHA256 and summarizes the full file. Neither mode filters reads.
"""

import argparse
import gzip
import hashlib
import json
import re
import statistics
import sys
import time
import urllib.request
import zlib
from pathlib import Path

REGISTRY = Path(__file__).with_name("long_read_runs.json")
CHUNK = 1 << 20

# PacBio read names: movie/zmw/start_end (subread), movie/zmw/ccs, movie/zmw.
SUBREAD = re.compile(r"^(m[^/\s]+)/(\d+)/(\d+)_(\d+)$")
CCS = re.compile(r"^(m[^/\s]+)/(\d+)/ccs(?:/.*)?$")
ZMW_READ = re.compile(r"^(m[^/\s]+)/(\d+)$")


def classify_name(header):
    """Return (read_unit, zmw_key, sra_spot) for an ENA/SRA FASTQ header line."""
    fields = header[1:].split()
    spot = fields[0].rsplit(".", 1)[-1] if fields else ""
    name = fields[1] if len(fields) > 1 else ""
    for unit, pattern in (("subread", SUBREAD), ("ccs", CCS), ("zmw_read", ZMW_READ)):
        match = pattern.match(name)
        if match:
            return unit, (match.group(1), int(match.group(2))), spot
    if not name or re.fullmatch(r"\d*/\d+", name):
        return "name_stripped", None, spot
    return "other_name", None, spot


class Summary:
    def __init__(self, platform):
        self.platform = platform
        self.spots = {}
        self.lengths = []
        self.units = {}
        self.zmws = set()
        self.placeholder_quality = 0
        self.examples = []

    def add(self, header, seq, qual):
        unit, zmw, spot = classify_name(header)
        self.spots[spot] = self.spots.get(spot, 0) + 1
        self.units[unit] = self.units.get(unit, 0) + 1
        if zmw is not None:
            self.zmws.add(zmw)
        if qual and qual.count(b"!") == len(qual):
            self.placeholder_quality += 1
        self.lengths.append(len(seq))
        if len(self.examples) < 4:
            self.examples.append(header)

    def report(self):
        n = len(self.lengths)
        total = sum(self.lengths)
        result = dict(
            records=n,
            bases=total,
            read_unit_counts=self.units,
            distinct_zmws=len(self.zmws) if self.zmws else None,
            distinct_sra_spots=len(self.spots),
            max_records_per_spot=max(self.spots.values(), default=0),
            placeholder_quality_records=self.placeholder_quality,
            example_headers=self.examples,
        )
        if n:
            ordered = sorted(self.lengths, reverse=True)
            acc = 0
            for length in ordered:
                acc += length
                if acc * 2 >= total:
                    break
            result.update(
                mean_length=round(total / n, 1),
                median_length=statistics.median(ordered),
                max_length=ordered[0],
                n50=length,
            )
        result["molecule_grouping"] = self.grouping(n)
        return result

    def grouping(self, n):
        """How records map to sequenced molecules; never inferred from length alone."""
        if self.platform == "OXFORD_NANOPORE":
            return "one_record_per_molecule"
        if n and self.units.get("subread", 0) + self.units.get("ccs", 0) + self.units.get("zmw_read", 0) == n:
            return "zmw_from_read_name"
        if max(self.spots.values(), default=0) > 1:
            # SRA loads PacBio RS runs with one spot per ZMW; confirm against SRA metadata.
            return "sra_spot_as_zmw_assumed"
        return "molecule_grouping_unavailable"


def parse_fastq(lines, summary):
    """Feed complete 4-line records; return the number consumed."""
    records = len(lines) // 4
    for i in range(records):
        header, seq, plus, qual = lines[4 * i : 4 * i + 4]
        if not header.startswith(b"@") or not plus.startswith(b"+"):
            raise ValueError("Malformed FASTQ record near: " + header[:80].decode(errors="replace"))
        summary.add(header.decode(errors="replace"), seq, qual)
    return records


def probe(run, nbytes):
    request = urllib.request.Request(run["url"], headers={"Range": "bytes=0-%d" % (nbytes - 1)})
    with urllib.request.urlopen(request, timeout=120) as response:
        data = response.read()
    text = zlib.decompressobj(zlib.MAX_WBITS | 16).decompress(data)
    lines = text.split(b"\n")[:-1]  # last line may be truncated
    summary = Summary(run["platform"])
    parse_fastq(lines[: len(lines) // 4 * 4], summary)
    result = summary.report()
    result["compressed_bytes_read"] = len(data)
    result["scope"] = "probe_prefix_only"
    return result


def download(url, target, expected_bytes):
    part = target.with_name(target.name + ".part")
    offset = part.stat().st_size if part.exists() else 0
    if offset > expected_bytes:
        raise SystemExit("Partial file larger than archive size: " + str(part))
    if offset < expected_bytes:
        headers = {"Range": "bytes=%d-" % offset} if offset else {}
        with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=300) as response:
            if offset and response.status != 206:
                raise SystemExit("Server ignored resume request for " + url)
            with open(part, "ab") as handle:
                while True:
                    block = response.read(CHUNK)
                    if not block:
                        break
                    handle.write(block)
    part.rename(target)


def file_hashes(path):
    md5, sha256 = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def summarize_file(path, platform):
    summary = Summary(platform)
    with gzip.open(path, "rb") as handle:
        buffer = []
        for line in handle:
            buffer.append(line.rstrip(b"\r\n"))
            if len(buffer) == 4:
                parse_fastq(buffer, summary)
                buffer = []
        if buffer:
            raise ValueError("Truncated FASTQ: %d trailing lines" % len(buffer))
    return summary.report()


def fetch(run, out):
    work = out / run["run"]
    work.mkdir(parents=True, exist_ok=True)
    report_path = work / "verify.json"
    if report_path.exists():
        raise SystemExit("Refusing to overwrite existing " + str(report_path))
    target = work / run["url"].rsplit("/", 1)[1]
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    if not target.exists():
        download(run["url"], target, run["fastq_bytes"])
    size = target.stat().st_size
    md5, sha256 = file_hashes(target)
    checks = dict(
        bytes_match=size == run["fastq_bytes"],
        md5_match=md5 == run["fastq_md5"],
    )
    result = dict(run=run, download_started=started, local_file=str(target), bytes=size, md5=md5, sha256=sha256)
    if all(checks.values()):
        stats = summarize_file(target, run["platform"])
        # ENA read_count counts SRA spots; for PacBio RS a spot can hold several subreads.
        checks.update(
            spot_count_match=stats["distinct_sra_spots"] == run["archive_read_count"],
            base_count_match=stats["bases"] == run["archive_base_count"],
        )
        result["summary"] = stats
    result["checks"] = checks
    result["status"] = "verified" if all(checks.values()) else "mismatch"
    report_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["probe", "fetch"])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--runs", nargs="*", help="Run or sample IDs; default all registered runs")
    parser.add_argument("--probe-bytes", type=int, default=4 << 20)
    args = parser.parse_args()
    runs = json.loads(REGISTRY.read_text())["runs"]
    if args.runs:
        wanted = set(args.runs)
        runs = [r for r in runs if r["run"] in wanted or r["sample_id"] in wanted]
        missing = wanted - {r["run"] for r in runs} - {r["sample_id"] for r in runs}
        if missing:
            raise SystemExit("Unknown runs: " + ", ".join(sorted(missing)))
    args.out.mkdir(parents=True, exist_ok=True)
    failed = False
    for run in runs:
        if args.mode == "probe":
            stats = probe(run, args.probe_bytes)
            (args.out / (run["run"] + ".probe.json")).write_text(json.dumps(stats, indent=2) + "\n")
            status = "probed"
        else:
            result = fetch(run, args.out)
            stats = result.get("summary", {})
            status = result["status"]
            failed |= status != "verified"
        print(
            run["sample_id"],
            run["run"],
            status,
            "records=%s" % stats.get("records"),
            "units=%s" % stats.get("read_unit_counts"),
            "zmws=%s" % stats.get("distinct_zmws"),
            "spots=%s" % stats.get("distinct_sra_spots"),
            "max_per_spot=%s" % stats.get("max_records_per_spot"),
            "placeholder_qual=%s" % stats.get("placeholder_quality_records"),
            stats.get("molecule_grouping"),
        )
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
