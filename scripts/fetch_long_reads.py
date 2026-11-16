"""Download and verify the registered long-read runs (docs/LONG_READ_AVAILABILITY.md).

`probe` reads only the first bytes of each archive FASTQ to classify read names,
read unit and quality encoding. `fetch` downloads, checks archive size/MD5,
records SHA256 and summarizes the full file. Neither mode filters reads.

Both modes use the network and refuse to run inside a SLURM job: downloads belong
to the transfer node (scripts/reads_ready.py prepare), compute jobs only check
(scripts/reads_ready.py check). A download stays <file>.part while it transfers
and <file>.unverified while it is checked; only a run whose checks all pass gets
its final file name and verify.json status verified.
"""

import argparse
import gzip
import hashlib
import json
import re
import statistics
import subprocess
import sys
import time
import urllib.request
import zlib
from pathlib import Path

from reads_ready import refuse_in_slurm

REGISTRY = Path(__file__).with_name("long_read_runs.json")
CHUNK = 1 << 20
MD5_ATTEMPTS = 3

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
    refuse_in_slurm("probe")
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


def download(url, target, expected_bytes, attempts=8):
    # curl, not urllib: urllib downloads of full length repeatedly carried a
    # short corrupted stretch, while curl copies of the same files were clean.
    refuse_in_slurm()
    part = target.with_name(target.name + ".part")
    for attempt in range(attempts):
        offset = part.stat().st_size if part.exists() else 0
        if offset > expected_bytes:
            raise SystemExit("Partial file larger than archive size: " + str(part))
        if offset == expected_bytes:
            break
        cmd = ["curl", "-sS", "-f", "-L", "--connect-timeout", "60", "--speed-limit", "1000",
               "--speed-time", "300", "-C", "-", "-o", str(part), url]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("download interrupted (curl exit %d: %s); attempt %d/%d"
                  % (result.returncode, result.stderr.strip(), attempt + 1, attempts), file=sys.stderr)
            time.sleep(30)
    if not part.exists() or part.stat().st_size != expected_bytes:
        raise SystemExit("Download incomplete after %d attempts: %s" % (attempts, part))
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


def fetch(run, out, content=True):
    """Download one run and check it. content=False stops after bytes and MD5 (status md5_verified), for a
    transfer node; check_content then completes it (no network). Only a full check gives status verified."""
    work = out / run["run"]
    work.mkdir(parents=True, exist_ok=True)
    report_path = work / "verify.json"
    if report_path.exists():
        raise SystemExit("Refusing to overwrite existing " + str(report_path))
    target = work / run["url"].rsplit("/", 1)[1]
    pending = unverified(target)
    if target.exists() and not pending.exists():
        target.rename(pending)  # final name without verify.json (interrupted older fetch): check it again
    old_part, part = target.with_name(target.name + ".part"), pending.with_name(pending.name + ".part")
    if old_part.exists() and not part.exists():
        old_part.rename(part)  # resume a partial download made before the .unverified naming
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    # Full-length files with a short corrupted stretch have been seen; keep each
    # bad copy for the record and download again.
    rejected = []
    for attempt in range(MD5_ATTEMPTS):
        if not pending.exists():
            download(run["url"], pending, run["fastq_bytes"])
        size = pending.stat().st_size
        md5, sha256 = file_hashes(pending)
        if md5 == run["fastq_md5"] or attempt == MD5_ATTEMPTS - 1:
            break
        bad = work / ("md5_mismatch_attempt%d" % (attempt + 1))
        bad.mkdir(exist_ok=True)
        pending.rename(bad / target.name)
        rejected.append(dict(attempt=attempt + 1, md5=md5, bytes=size, kept_as=str(bad / target.name)))
    checks = dict(
        bytes_match=size == run["fastq_bytes"],
        md5_match=md5 == run["fastq_md5"],
    )
    result = dict(run=run, download_started=started, local_file=str(pending), bytes=size, md5=md5, sha256=sha256,
                  rejected_downloads=rejected)
    if all(checks.values()) and content:
        content_checks(result, pending, run, checks)
    result["checks"] = checks
    if not all(checks.values()):
        result["status"] = "mismatch"
    else:
        result["status"] = "verified" if content else "md5_verified"
    if result["status"] == "verified":
        pending.rename(target)
        result["local_file"] = str(target)
    report_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def unverified(target):
    return target.with_name(target.name + ".unverified")


def content_checks(result, target, run, checks):
    stats = summarize_file(target, run["platform"])
    # ENA read_count counts SRA spots; for PacBio RS a spot can hold several subreads.
    # Spots absent from the FASTQ are accepted only if they carry no bases,
    # i.e. the base count still matches exactly (empty ZMWs are not exported).
    bases_match = stats["bases"] == run["archive_base_count"]
    missing = run["archive_read_count"] - stats["distinct_sra_spots"]
    checks.update(
        spot_count_match=missing == 0 or (missing > 0 and bases_match),
        base_count_match=bases_match,
    )
    result["spots_absent_from_fastq"] = missing
    result["summary"] = stats


def check_content(work):
    """Complete an md5_verified download in place: the file must still hash as recorded, then the FASTQ
    content checks run as in fetch. Any other status is returned unchanged."""
    report_path = work / "verify.json"
    result = json.loads(report_path.read_text())
    if result["status"] != "md5_verified":
        return result
    target = Path(result["local_file"])
    checks = result["checks"]
    final = target.with_name(target.name[: -len(".unverified")]) if target.name.endswith(".unverified") else target
    md5, sha256 = file_hashes(target)
    checks["unchanged_since_download"] = md5 == result["md5"] and sha256 == result["sha256"]
    if all(checks.values()):
        content_checks(result, target, result["run"], checks)
    result["content_checked"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    result["status"] = "verified" if all(checks.values()) else "mismatch"
    if result["status"] == "verified" and final != target:
        target.rename(final)
        result["local_file"] = str(final)
    report_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def select_runs(ids):
    """Registry runs by run or sample id (all registered runs when ids is empty)."""
    runs = json.loads(REGISTRY.read_text())["runs"]
    if ids:
        wanted = set(ids)
        runs = [r for r in runs if r["run"] in wanted or r["sample_id"] in wanted]
        missing = wanted - {r["run"] for r in runs} - {r["sample_id"] for r in runs}
        if missing:
            raise SystemExit("Unknown runs: " + ", ".join(sorted(missing)))
    return runs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["probe", "fetch"])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--runs", nargs="*", help="Run or sample IDs; default all registered runs")
    parser.add_argument("--probe-bytes", type=int, default=4 << 20)
    args = parser.parse_args()
    refuse_in_slurm(args.mode)
    runs = select_runs(args.runs)
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
