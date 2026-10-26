"""Pilot of other ENA transports (Kingfisher ena-ascp / ena-ftp) against the curl fetch in fetch_long_reads.py.

Each run is downloaded into <out>/<method>/<run>/ (never into the shared reads directory) and checked exactly as
fetch_long_reads.fetch does: archive bytes and MD5, then the whole FASTQ summarized for base and spot counts.
report.json per run and method keeps tool versions, the command, wall time, return code and the checks; the
Kingfisher log is kept next to it. Run metadata comes from the pool's runs.tsv, as for the production fetch.
Usage: python ena_transport_pilot.py --pool-dir <pool> --out <dir> --runs R1 [R2 R3] --methods ena-ascp [ena-ftp]
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
import fetch_long_reads as flr  # noqa: E402
from s3_readcheck import read  # noqa: E402


def version(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return (r.stdout + r.stderr).strip().splitlines()[0] if (r.stdout + r.stderr).strip() else "rc=%d" % r.returncode
    except Exception as e:
        return "unavailable: %r" % e


def one(r, method, out, key=None):
    work = out / method / r["run_accession"]
    if (work / "report.json").exists():
        raise SystemExit("Refusing to overwrite existing " + str(work / "report.json"))
    work.mkdir(parents=True, exist_ok=True)
    cmd = ["kingfisher", "get", "-r", r["run_accession"], "-m", method, "-f", "fastq.gz", "--output-directory", str(work)]
    if key and method == "ena-ascp":
        cmd += ["--ascp-ssh-key", str(key)]
    started = time.time()
    with open(str(work / "kingfisher.log"), "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.time() - started
    expected = r["fastq_ftp"].rsplit("/", 1)[1]
    target = work / expected
    found = sorted(p.name for p in work.iterdir() if p.name not in ("kingfisher.log", "report.json"))
    run = dict(run=r["run_accession"], fastq_bytes=int(r["fastq_bytes"]), fastq_md5=r["fastq_md5"],
               platform=r["instrument_platform"], archive_base_count=int(r["base_count"]), archive_read_count=int(r["read_count"]))
    result = dict(run=run, method=method, command=cmd, return_code=rc, wall_seconds=round(wall, 1), files=found,
                  expected_file=expected)
    checks = dict(file_present=target.exists())
    if target.exists():
        size = target.stat().st_size
        md5, sha256 = flr.file_hashes(target)
        checks.update(bytes_match=size == run["fastq_bytes"], md5_match=md5 == run["fastq_md5"])
        result.update(bytes=size, md5=md5, sha256=sha256, mb_per_s=round(size / 1e6 / wall, 2) if wall else None)
        if checks["bytes_match"] and checks["md5_match"]:
            stats = flr.summarize_file(target, run["platform"])
            bases_match = stats["bases"] == run["archive_base_count"]
            missing = run["archive_read_count"] - stats["distinct_sra_spots"]
            checks.update(spot_count_match=missing == 0 or (missing > 0 and bases_match), base_count_match=bases_match)
            result["spots_absent_from_fastq"] = missing
            result["summary"] = stats
    result["checks"] = checks
    result["status"] = "verified" if all(checks.values()) and rc == 0 else "failed"
    (work / "report.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pool-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--runs", nargs="+", required=True)
    p.add_argument("--ascp-ssh-key", type=Path, help="key for ascp -i (default: the one bundled with Kingfisher)")
    p.add_argument("--methods", nargs="+", default=["ena-ascp"], choices=["ena-ascp", "ena-ftp"])
    a = p.parse_args()
    if len(a.runs) > 3:
        raise SystemExit("pilot: at most three runs")
    runs = {r["run_accession"]: r for r in read(a.pool_dir / "runs.tsv")}
    missing = [x for x in a.runs if x not in runs]
    if missing:
        raise SystemExit("Unknown runs: " + ", ".join(missing))
    a.out.mkdir(parents=True, exist_ok=True)
    env = dict(host=subprocess.run(["hostname"], capture_output=True, text=True).stdout.strip(),
               kingfisher=version(["kingfisher", "--version"]), ascp=version(["ascp", "--version"]),
               curl=version(["curl", "--version"]), ascp_path=shutil.which("ascp"),
               ascp_ssh_key=str(a.ascp_ssh_key) if a.ascp_ssh_key else "kingfisher bundled")
    (a.out / ("tools_%s.json" % time.strftime("%Y%m%dT%H%M%S"))).write_text(json.dumps(env, indent=1) + "\n")
    print(json.dumps(env))
    for x in a.runs:
        for m in a.methods:
            res = one(runs[x], m, a.out, a.ascp_ssh_key)
            print(x, m, res["status"], "rc=%d" % res["return_code"], "%ss" % res["wall_seconds"], res.get("mb_per_s"), "MB/s",
                  json.dumps(res["checks"]), flush=True)


if __name__ == "__main__":
    main()
