"""Long-read inputs: downloaded and verified on the transfer node, only checked on compute nodes.

prepare  Transfer node (lrc-xfer, outside SLURM; refuses to run inside a SLURM job). Each run is downloaded with
         resume into <file>.part, kept as <file>.unverified while size, MD5 and FASTQ content are checked
         (fetch_long_reads.fetch), and gets its final name and verify.json status verified only when every check
         passes. Runs already verified are kept as they are; a run left md5_verified by an earlier transfer-only
         download gets its content check. Runs come from the registry (--runs: run or sample ids) or from an ENA run
         table (--runs-tsv: run_accession, sample_accession, fastq_ftp, fastq_bytes, fastq_md5,
         instrument_platform, base_count, read_count; --only limits it to some runs). Ends by rewriting the ready
         manifest. Storage is checked per batch: before anything is downloaded, prepare prints the reads directory
         usage, the bytes this batch still needs (runs not yet verified, minus any .part/.unverified already on
         disk) and the free space of the filesystem, and refuses the batch (exit 4) if less than --min-free-gb would
         remain afterwards (--dry-run stops after this check). Usage is printed again at the end.
index    Rewrite the ready manifest from the verify.json files under --reads (no network).
check    Compute side. Every requested run must be in the ready manifest with status verified, its file must exist
         with the recorded size (and SHA256 with --rehash), and its verify.json must still be verified with that
         SHA256. Otherwise one `input_not_ready <run> <reason>` line per run and exit 3. Never downloads; this module
         imports the downloader only inside `prepare`.

Ready manifest (<reads>/ready_manifest.tsv): one row per run directory with a verify.json (run, sample,
local_file, bytes, md5, sha256, status, verified_at, verify_json); only rows with status verified are usable.
"""

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

READS = Path("/global/scratch/users/kh36969/longread_callout/data/reads")
MANIFEST = "ready_manifest.tsv"
COLS = ["run", "sample", "local_file", "bytes", "md5", "sha256", "status", "verified_at", "verify_json"]
NOT_READY = 3
NO_SPACE = 4


class InputNotReady(SystemExit):
    """Raised (exit code 3) when a compute step is asked to use runs that are not ready."""

    def __init__(self, reasons):
        self.reasons = reasons
        lines = "\n".join("input_not_ready %s %s" % kv for kv in sorted(reasons.items()))
        print(lines, file=sys.stderr, flush=True)
        super().__init__(NOT_READY)


def in_slurm():
    return bool(os.environ.get("SLURM_JOB_ID"))


def refuse_in_slurm(what="download"):
    if in_slurm():
        raise SystemExit("%s refused inside SLURM job %s: downloads run on the transfer node (lrc-xfer) with "
                         "scripts/reads_ready.py prepare; compute jobs only read verified local files"
                         % (what, os.environ["SLURM_JOB_ID"]))


def index(reads=READS):
    rows = []
    for v in sorted(Path(reads).glob("*/verify.json")):
        j = json.loads(v.read_text())
        run = j.get("run", {})
        rows.append(dict(run=v.parent.name, sample=run.get("sample_id", ""), local_file=j.get("local_file", ""),
                         bytes=j.get("bytes", ""), md5=j.get("md5", ""), sha256=j.get("sha256", ""), status=j["status"],
                         verified_at=j.get("content_checked") or j.get("download_started", ""), verify_json=str(v)))
    out = Path(reads) / MANIFEST
    tmp = out.with_name(out.name + ".tmp")
    with open(str(tmp), "w") as h:
        w = csv.DictWriter(h, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.rename(out)
    return rows


def load(manifest):
    with open(str(manifest)) as h:
        return {r["run"]: r for r in csv.DictReader(h, delimiter="\t")}


def sha256_of(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def not_ready(runs, reads=READS, manifest=None, rehash=False):
    """{run: reason} for every requested run that may not be used; empty when all are ready."""
    path = Path(manifest) if manifest else Path(reads) / MANIFEST
    if not path.exists():
        return {r: "ready_manifest_missing:%s" % path for r in runs}
    m = load(path)
    bad = {}
    for run in runs:
        r = m.get(run)
        if r is None:
            bad[run] = "not_in_ready_manifest"
            continue
        if r["status"] != "verified":
            bad[run] = "status_%s" % r["status"]
            continue
        f = Path(r["local_file"])
        if not f.is_file():
            bad[run] = "file_missing:%s" % f
            continue
        if str(f.stat().st_size) != r["bytes"]:
            bad[run] = "size_mismatch:%d!=%s" % (f.stat().st_size, r["bytes"])
            continue
        v = Path(r["verify_json"])
        if not v.is_file():
            bad[run] = "verify_json_missing"
            continue
        j = json.loads(v.read_text())
        if j.get("status") != "verified" or j.get("sha256") != r["sha256"] or j.get("local_file") != r["local_file"]:
            bad[run] = "verify_json_differs_from_manifest"
            continue
        if rehash and sha256_of(f) != r["sha256"]:
            bad[run] = "sha256_mismatch"
    return bad


def require(runs, reads=READS, manifest=None, rehash=False):
    """{run: local_file} of runs that are all ready; raises InputNotReady (exit 3) otherwise."""
    runs = list(runs)
    bad = not_ready(runs, reads, manifest, rehash)
    if bad:
        raise InputNotReady(bad)
    path = Path(manifest) if manifest else Path(reads) / MANIFEST
    m = load(path)
    return {r: m[r]["local_file"] for r in runs}


def storage(reads):
    """(bytes of all files, of verified files, of .part/.unverified files, free bytes on the filesystem) of the
    reads dir; all files include rejected copies (md5_mismatch_attempt*) and logs."""
    reads = Path(reads)
    total = sum(f.stat().st_size for f in reads.rglob("*") if f.is_file())
    done = partial = 0
    for v in reads.glob("*/verify.json"):
        f = Path(json.loads(v.read_text()).get("local_file", ""))
        done += f.stat().st_size if f.is_file() else 0
    for f in list(reads.glob("*/*.part")) + list(reads.glob("*/*.unverified")):
        partial += f.stat().st_size
    s = os.statvfs(str(reads))
    return total, done, partial, s.f_bavail * s.f_frsize


def report_storage(reads, label):
    total, done, partial, free = storage(reads)
    print("storage %s: reads dir %.1f GB in all (%.1f GB verified, %.1f GB partial); filesystem free %.1f TB"
          % (label, total / 1e9, done / 1e9, partial / 1e9, free / 1e12), flush=True)
    return total, done, partial, free


def check_space(runs, reads, min_free_gb):
    """Refuse (exit 4) when this batch would leave less than min_free_gb free."""
    need = 0
    for run in runs:
        d = Path(reads) / run["run"]
        if (d / "verify.json").exists():
            continue
        have = sum(f.stat().st_size for f in list(d.glob("*.part")) + list(d.glob("*.unverified")))
        need += max(int(run["fastq_bytes"]) - have, 0)
    free = report_storage(reads, "before")[-1]
    print("batch: %d runs, %.2f GB still to download; free after %.1f GB (minimum %.0f GB)"
          % (len(runs), need / 1e9, (free - need) / 1e9, min_free_gb), flush=True)
    if free - need < min_free_gb * 1e9:
        print("prepare refused: not enough free space for this batch", file=sys.stderr)
        sys.exit(NO_SPACE)


def cmd_prepare(a):
    refuse_in_slurm("prepare")
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fetch_long_reads as flr  # downloader: imported here only, never by the check path
    if a.runs_tsv:
        with open(str(a.runs_tsv)) as h:
            rows = list(csv.DictReader(h, delimiter="\t"))
        only = set(a.only or [])
        runs = [dict(run=r["run_accession"], sample_id=r["sample_accession"], url="https://" + r["fastq_ftp"],
                     fastq_bytes=int(r["fastq_bytes"]), fastq_md5=r["fastq_md5"], platform=r["instrument_platform"],
                     archive_base_count=int(r["base_count"]), archive_read_count=int(r["read_count"]))
                for r in rows if not only or r["run_accession"] in only]
        multi = [r["run"] for r in runs if ";" in r["url"]]
        for r in multi:
            print(r, "skipped: more than one FASTQ file", file=sys.stderr)
        runs = [r for r in runs if r["run"] not in multi]
    else:
        runs = flr.select_runs(a.runs)
    a.reads.mkdir(parents=True, exist_ok=True)
    check_space(runs, a.reads, a.min_free_gb)
    if a.dry_run:
        return
    bad = []
    for run in runs:
        work = a.reads / run["run"]
        v = work / "verify.json"
        try:
            if v.exists():
                st = json.loads(v.read_text())["status"]
                if st == "md5_verified":
                    st = flr.check_content(work)["status"]
            else:
                st = flr.fetch(run, a.reads)["status"]
        except (Exception, SystemExit) as e:  # reported per run; a rerun resumes the .part file
            st = "error: %r" % e
        print(run["sample_id"], run["run"], st, flush=True)
        if st != "verified":
            bad.append(run["run"])
    rows = index(a.reads)
    report_storage(a.reads, "after")
    print("ready manifest: %d runs, %d verified; %d of %d requested runs not verified%s" % (
        len(rows), sum(r["status"] == "verified" for r in rows), len(bad), len(runs),
        (": " + " ".join(bad)) if bad else ""))
    sys.exit(1 if bad else 0)


def cmd_index(a):
    rows = index(a.reads)
    print("ready manifest %s: %d runs, %d verified" % (a.reads / MANIFEST, len(rows),
                                                       sum(r["status"] == "verified" for r in rows)))


def cmd_check(a):
    runs = list(a.runs or [])
    if a.runs_file:
        with open(str(a.runs_file)) as h:
            rows = list(csv.DictReader(h, delimiter="\t"))
        runs += [r[a.run_column] for r in rows if not (a.skip_multi_file and ";" in r.get("fastq_ftp", ""))]
    if not runs:
        raise SystemExit("check: no runs given")
    files = require(runs, a.reads, a.manifest, a.rehash)
    for r in runs:
        print("ready", r, files[r])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["prepare", "index", "check"])
    p.add_argument("--reads", type=Path, default=READS)
    p.add_argument("--runs", nargs="*", help="prepare: registry run or sample ids; check: run ids")
    p.add_argument("--runs-tsv", type=Path, help="prepare: ENA run table instead of the registry")
    p.add_argument("--only", nargs="*", help="prepare with --runs-tsv: only these run accessions")
    p.add_argument("--min-free-gb", type=float, default=500,
                   help="prepare: refuse the batch if less than this would stay free on the filesystem")
    p.add_argument("--dry-run", action="store_true", help="prepare: only the storage check, no download")
    p.add_argument("--runs-file", type=Path, help="check: TSV with a run column")
    p.add_argument("--run-column", default="run")
    p.add_argument("--skip-multi-file", action="store_true",
                   help="check with --runs-file: leave out runs whose fastq_ftp lists more than one file")
    p.add_argument("--manifest", type=Path, help="check: ready manifest (default <reads>/%s)" % MANIFEST)
    p.add_argument("--rehash", action="store_true", help="check: recompute SHA256 of every file")
    a = p.parse_args()
    dict(prepare=cmd_prepare, index=cmd_index, check=cmd_check)[a.mode](a)


if __name__ == "__main__":
    main()
