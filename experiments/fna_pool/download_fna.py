"""S1: resumable NCBI Datasets download and verification of the FNA/GFF3/FAA manifest written by build_pool.py.

run     Takes manifest rows not yet verified (failed rows are retried up to --max-attempts), in chunks:
        dehydrated package -> rehydrate -> per-accession checks -> move into <data>/genomes/<accession>/.
        A file only becomes usable after its accession passes; every attempt is appended to <data>/state.jsonl.
        Re-running continues where it stopped. One instance at a time (<data>/.lock).
report  Pair check (GCF vs GCA record sequences) and summary: <data>/report.json, plus small tables for git
        (--out): files.tsv (bytes, MD5, SHA256 per file), accessions.tsv (checks per accession), pair_checks.tsv.

Checks per accession: package MD5 for every file; FNA records equal the NCBI sequence report (accession, length),
all records kept; assembly report accession and BioSample equal the manifest; GFF3 sequence-region ids/lengths
equal the FNA; every GFF3 CDS protein_id is in the FAA and every FAA id is in the GFF3. No annotation in the
package is recorded as annotation_missing, not as a failure. Failures are never counted as zero hits.

Python 3.6+, standard library only, plus the datasets CLI. Meant to run on lrc-xfer (see README.md).
"""

import argparse
import csv
import glob
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
import time
import zipfile
from collections import Counter, OrderedDict
from pathlib import Path

DATASETS = "/global/home/users/kh36969/.conda/envs/claude-env/bin/datasets"
PART = {"genome": "_genomic.fna", "gff3": "genomic.gff", "protein": "protein.faa", "seq-report": "sequence_report.jsonl"}


def log(*a):
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), *a, flush=True)


def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(str(path), "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            md5.update(b)
            sha.update(b)
    return md5.hexdigest(), sha.hexdigest()


def load_state(data):
    state = {}
    p = data / "state.jsonl"
    if p.exists():
        for line in open(str(p)):
            e = json.loads(line)
            s = state.setdefault(e["accession"], dict(attempts=0))
            s["status"] = e["status"]
            if e["status"] in ("verified", "failed"):
                s["attempts"] += 1
            s["last"] = e
    return state


def append_state(data, **e):
    e["time"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    with open(str(data / "state.jsonl"), "a") as h:
        h.write(json.dumps(e) + "\n")


def read_fasta(path):
    """id -> (length, SHA256 of upper-case sequence)."""
    out, name, h, n = OrderedDict(), None, None, 0
    with open(str(path), "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if name is not None:
                    out[name] = (n, h.hexdigest())
                name, h, n = line[1:].split()[0].decode(), hashlib.sha256(), 0
            else:
                s = line.strip().upper()
                h.update(s)
                n += len(s)
    if name is not None:
        out[name] = (n, h.hexdigest())
    return out


def read_gff(path):
    regions, proteins, cds, pseudo, partial = {}, set(), 0, 0, 0
    with open(str(path)) as f:
        for line in f:
            if line.startswith("##sequence-region"):
                _, sid, start, end = line.split()[:4]
                regions[sid] = int(end) - int(start) + 1
            elif line.startswith("#") or not line.strip():
                continue
            else:
                col = line.rstrip("\n").split("\t")
                if len(col) < 9 or col[2] != "CDS":
                    continue
                cds += 1
                attrs = dict(kv.split("=", 1) for kv in col[8].split(";") if "=" in kv)
                if "protein_id" in attrs:
                    proteins.add(attrs["protein_id"])
                pseudo += attrs.get("pseudo") == "true"
                partial += attrs.get("partial") == "true"
    return dict(regions=regions, proteins=proteins, cds_rows=cds, pseudo_cds_rows=pseudo, partial_cds_rows=partial)


def faa_ids(path):
    with open(str(path)) as f:
        return {line[1:].split()[0] for line in f if line.startswith(">")}


def verify(acc, row, pkg, md5ok):
    d = pkg / "ncbi_dataset" / "data" / acc
    if not d.is_dir():
        return dict(status="failed", reason="not_in_package"), []
    want = row["include"].split(",")
    files = sorted(p for p in d.iterdir() if p.is_file())
    checks, info = OrderedDict(), OrderedDict()
    checks["package_md5_all_ok"] = all(md5ok.get("ncbi_dataset/data/%s/%s" % (acc, p.name), p.name == "sequence_report.jsonl")
                                       for p in files)
    rep = None
    for line in open(str(pkg / "ncbi_dataset" / "data" / "assembly_data_report.jsonl")):
        r = json.loads(line)
        if r.get("accession") == acc:
            rep = r
    checks["assembly_report_present"] = rep is not None
    if rep:
        info.update(biosample=rep["assemblyInfo"]["biosample"]["accession"], assembly_level=rep["assemblyInfo"]["assemblyLevel"],
                    total_bp=int(rep["assemblyStats"]["totalSequenceLength"]),
                    annotation=rep.get("annotationInfo", {}).get("name", ""),
                    annotation_date=rep.get("annotationInfo", {}).get("releaseDate", ""))
        checks["biosample_matches_manifest"] = info["biosample"] == row["biosample"]
    fna = [p for p in files if p.name.endswith(PART["genome"])]
    checks["one_fna"] = len(fna) == 1
    if len(fna) == 1:
        recs = read_fasta(fna[0])
        key = "refseqAccession" if acc.startswith("GCF_") else "genbankAccession"
        sr = [json.loads(l) for l in open(str(d / PART["seq-report"]))] if (d / PART["seq-report"]).exists() else []
        # Some NCBI sequence-report rows omit length; those records are matched by accession only
        # (total length and GFF3 sequence-region lengths are still checked below).
        srmap = {s.get(key): int(s["length"]) if "length" in s else None for s in sr}
        checks["seq_report_present"] = bool(sr)
        checks["fna_records_equal_seq_report"] = set(recs) == set(srmap) and all(
            srmap[k] is None or recs[k][0] == srmap[k] for k in recs)
        info["seq_report_rows_without_length"] = sum(1 for v in srmap.values() if v is None)
        checks["fna_total_equals_assembly_report"] = bool(rep) and sum(v[0] for v in recs.values()) == info["total_bp"]
        info.update(records=len(recs), fna_bp=sum(v[0] for v in recs.values()),
                    record_sha256={k: v[1] for k, v in recs.items()},
                    molecule_types=dict(Counter(s.get("assignedMoleculeLocationType", "") for s in sr)))
        if "gff3" in want:
            gff, faa = d / PART["gff3"], d / PART["protein"]
            if gff.exists() and faa.exists():
                g = read_gff(gff)
                ids = faa_ids(faa)
                checks["gff_regions_equal_fna"] = g["regions"] == {k: v[0] for k, v in recs.items()}
                checks["gff_protein_ids_in_faa"] = g["proteins"] <= ids
                checks["faa_ids_in_gff"] = ids <= g["proteins"]
                info.update(annotation_status="present", cds_rows=g["cds_rows"], pseudo_cds_rows=g["pseudo_cds_rows"],
                            partial_cds_rows=g["partial_cds_rows"], gff_protein_ids=len(g["proteins"]), faa_proteins=len(ids),
                            gff_ids_missing_from_faa=len(g["proteins"] - ids), faa_ids_missing_from_gff=len(ids - g["proteins"]))
            elif not gff.exists() and not faa.exists() and not info.get("annotation"):
                info["annotation_status"] = "annotation_missing"
            else:
                checks["annotation_files_complete"] = False
    status = "verified" if all(checks.values()) else "failed"
    failed = [k for k, v in checks.items() if not v]
    return dict(status=status, reason=";".join(failed), checks=checks, **info), files


def run(a):
    data = a.data
    (data / "genomes").mkdir(parents=True, exist_ok=True)
    (data / "tmp").mkdir(exist_ok=True)
    lock = data / ".lock"
    try:
        fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        sys.exit("Another run holds %s (%s); remove it only if that run is dead." % (lock, lock.read_text().strip()))
    os.write(fd, ("%s pid %d %s" % (socket.gethostname(), os.getpid(), time.strftime("%Y-%m-%dT%H:%M:%S"))).encode())
    os.close(fd)
    try:
        rows = list(csv.DictReader(open(str(a.manifest)), delimiter="\t"))
        state = load_state(data)
        todo = [r for r in rows if state.get(r["accession"], {}).get("status") != "verified"
                and state.get(r["accession"], {}).get("attempts", 0) < a.max_attempts]
        log("manifest %d rows; verified %d; to do %d" % (len(rows), sum(1 for r in rows if state.get(r["accession"], {}).get("status") == "verified"), len(todo)))
        groups = OrderedDict()
        for r in todo:
            groups.setdefault(r["include"], []).append(r)
        n = 0
        for include, rs in groups.items():
            for i in range(0, len(rs), a.chunk):
                if a.limit and n >= a.limit:
                    break
                chunk = rs[i:i + a.chunk][: (a.limit - n) if a.limit else None]
                n += len(chunk)
                do_chunk(a, data, include, chunk)
    finally:
        lock.unlink()
    log("done")


def do_chunk(a, data, include, chunk):
    tag = "chunk_%s_%d" % (time.strftime("%Y%m%d%H%M%S"), os.getpid())
    tmp = data / "tmp" / tag
    tmp.mkdir()
    (tmp / "acc.txt").write_text("".join(r["accession"] + "\n" for r in chunk))
    for r in chunk:
        append_state(data, accession=r["accession"], status="downloading", chunk=tag, host=socket.gethostname())
    log("%s: %d accessions, include=%s" % (tag, len(chunk), include))
    try:
        subprocess.run([DATASETS, "download", "genome", "accession", "--inputfile", str(tmp / "acc.txt"), "--include", include,
                        "--dehydrated", "--filename", str(tmp / "pkg.zip"), "--no-progressbar"], check=True,
                       stderr=subprocess.PIPE, universal_newlines=True)
        with zipfile.ZipFile(str(tmp / "pkg.zip")) as z:
            z.extractall(str(tmp / "pkg"))
        for attempt in range(3):
            p = subprocess.run([DATASETS, "rehydrate", "--directory", str(tmp / "pkg"), "--max-workers", str(a.workers),
                                "--no-progressbar"], stderr=subprocess.PIPE, universal_newlines=True)
            if p.returncode == 0:
                break
            log("rehydrate exit %d (%s); retry %d" % (p.returncode, p.stderr.strip()[-200:], attempt + 1))
            time.sleep(30)
    except (subprocess.CalledProcessError, zipfile.BadZipFile) as e:
        msg = getattr(e, "stderr", "") or str(e)
        for r in chunk:
            append_state(data, accession=r["accession"], status="failed", chunk=tag, reason="package: " + msg.strip()[-300:])
        log("%s failed: %s" % (tag, msg.strip()[-300:]))
        if not a.keep_failed:
            shutil.rmtree(str(tmp))
        return
    md5ok = {}
    for line in open(str(tmp / "pkg" / "md5sum.txt")):
        md5, rel = line.split(None, 1)
        rel = rel.strip()
        path = tmp / "pkg" / rel
        md5ok[rel] = path.exists() and hashes(path)[0] == md5
    for r in chunk:
        acc = r["accession"]
        res, files = verify(acc, r, tmp / "pkg", md5ok)
        if res["status"] == "verified":
            dest = data / "genomes" / acc
            if dest.exists():
                shutil.rmtree(str(dest.with_name(acc + ".old")), ignore_errors=True)
                dest.rename(dest.with_name(acc + ".old"))
            dest.mkdir()
            entries = []
            for f in files:
                md5, sha = hashes(f)
                entries.append(dict(file=f.name, bytes=f.stat().st_size, md5=md5, sha256=sha))
                shutil.move(str(f), str(dest / f.name))
            rep = [json.loads(l) for l in open(str(tmp / "pkg" / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"))]
            (dest / "assembly_data_report.json").write_text(json.dumps([x for x in rep if x.get("accession") == acc][0], indent=1) + "\n")
            res.update(accession=acc, biosample_manifest=r["biosample"], purpose=r["purpose"], include=r["include"], files=entries,
                       datasets_version=subprocess.run([DATASETS, "--version"], stdout=subprocess.PIPE,
                                                       universal_newlines=True).stdout.strip(), chunk=tag)
            (dest / "verify.json").write_text(json.dumps(res, indent=1) + "\n")
            shutil.rmtree(str(dest.with_name(acc + ".old")), ignore_errors=True)
        append_state(data, accession=acc, status=res["status"], chunk=tag, reason=res.get("reason", ""))
    ok = sum(1 for r in chunk if load_state(data)[r["accession"]]["status"] == "verified")
    log("%s: %d/%d verified" % (tag, ok, len(chunk)))
    if ok == len(chunk) or not a.keep_failed:
        shutil.rmtree(str(tmp))


def report(a):
    data = a.data
    rows = list(csv.DictReader(open(str(a.manifest)), delimiter="\t"))
    state = load_state(data)
    ver = {}
    for r in rows:
        p = data / "genomes" / r["accession"] / "verify.json"
        if state.get(r["accession"], {}).get("status") == "verified" and p.exists():
            ver[r["accession"]] = json.loads(p.read_text())
    pairs = []
    by_bs = {}
    for r in rows:
        by_bs.setdefault(r["biosample"], {})[r["purpose"]] = r["accession"]
    for bs, d in sorted(by_bs.items()):
        if "gca_identity_check" not in d:
            continue
        gcf, gca = d["primary"], d["gca_identity_check"]
        if gcf in ver and gca in ver:
            same = sorted(ver[gcf]["record_sha256"].values()) == sorted(ver[gca]["record_sha256"].values())
            status = "identical_sequences" if same else "different_sequences"
        else:
            status = "not_checked_download_incomplete"
        pairs.append(dict(biosample=bs, refseq=gcf, genbank=gca, status=status))
    files, accs = [], []
    for r in rows:
        acc, s = r["accession"], state.get(r["accession"], {})
        v = ver.get(acc, {})
        accs.append(OrderedDict(
            accession=acc, biosample=r["biosample"], purpose=r["purpose"], status=s.get("status", "pending"),
            attempts=s.get("attempts", 0), reason=s.get("last", {}).get("reason", ""), records=v.get("records", ""),
            fna_bp=v.get("fna_bp", ""), molecule_types=";".join("%s:%d" % kv for kv in sorted(v.get("molecule_types", {}).items())),
            annotation_status=v.get("annotation_status", ""), annotation=v.get("annotation", ""),
            annotation_date=v.get("annotation_date", ""), cds_rows=v.get("cds_rows", ""), faa_proteins=v.get("faa_proteins", "")))
        for f in v.get("files", []):
            files.append(OrderedDict(accession=acc, **f))
    out = dict(
        manifest=str(a.manifest), rows=len(rows), status=dict(Counter(x["status"] for x in accs)),
        by_purpose={p: dict(Counter(x["status"] for x in accs if x["purpose"] == p)) for p in sorted({x["purpose"] for x in accs})},
        annotation=dict(Counter(x["annotation_status"] for x in accs if x["purpose"] == "primary" and x["status"] == "verified")),
        pair_checks=dict(Counter(p["status"] for p in pairs)), bytes=sum(f["bytes"] for f in files), files=len(files),
        failed=[OrderedDict((k, x[k]) for k in ("accession", "attempts", "reason")) for x in accs if x["status"] == "failed"])
    (data / "report.json").write_text(json.dumps(out, indent=2) + "\n")
    if a.out:
        a.out.mkdir(parents=True, exist_ok=True)
        for name, tab in (("accessions.tsv", accs), ("files.tsv", files), ("pair_checks.tsv", pairs)):
            if tab:
                with open(str(a.out / name), "w") as h:
                    w = csv.DictWriter(h, list(tab[0]), delimiter="\t", lineterminator="\n")
                    w.writeheader()
                    w.writerows(tab)
        (a.out / "s1_report.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "failed"}, indent=2), "\nfailed:", len(out["failed"]))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["run", "report"])
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--data", type=Path, required=True, help="data root on scratch")
    p.add_argument("--out", type=Path, help="report: repo directory for small tables")
    p.add_argument("--chunk", type=int, default=50)
    p.add_argument("--workers", type=int, default=8, help="rehydrate download workers (1-30)")
    p.add_argument("--max-attempts", type=int, default=3)
    p.add_argument("--limit", type=int, default=0, help="stop after this many accessions (pilot)")
    p.add_argument("--keep-failed", action="store_true", help="keep tmp chunk dirs that had failures")
    a = p.parse_args()
    run(a) if a.mode == "run" else report(a)


if __name__ == "__main__":
    main()
