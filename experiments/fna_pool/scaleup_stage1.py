"""Scale-up stage 1 (docs/PLAN_SCALE_UP_STAGE1_ASSEMBLY_DISCOVERY.md): assembly-first discovery over frozen taxon batches.

Reuses build_pool / download_fna / build_s2 / s3_pair_calls / s4_annotate; this file holds the scale-up entry points.

baseline  Snapshot before the run: checkout, uncommitted changes, queued jobs, data directories, tool versions and
          binary hashes, frozen parameters of S2/S3/S4/read check, and hashes of the Comamonadaceae baseline results.
          Writes <out>/baseline.json.
targets   Wave definition from the archived priority catalog (--priorities, priority P1_lower_attention_review in
          review_rank order) joined to the assembly_first_v1 pools built per family under <out>/<family>/. Writes
          target_batches.tsv (one row per logical batch; reads-metadata counts kept apart from NCBI assembly
          counts), cross_batch_overlap.tsv (BioSamples or accessions in more than one batch) and
          baseline_pool_delta.tsv (old Comamonadaceae reads_linked_v1 pool vs the new pool, with reasons).
"""

import argparse
import hashlib
import importlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PROJ = Path("/global/scratch/users/kh36969/longread_callout")
TOOLS = {
    "minimap2": PROJ / "env/msi/bin/minimap2",
    "samtools": PROJ / "env/msi/bin/samtools",
    "python_msi": PROJ / "env/msi/bin/python",
    "mmseqs": PROJ / "env/s4/bin/mmseqs",
    "skani": PROJ / "paraburkholderia_ani/skani",
    "datasets": Path("/global/home/users/kh36969/.conda/envs/claude-env/bin/datasets"),
    "python_system": Path("/usr/bin/python3"),
}
VERSION_ARGS = {"mmseqs": ["version"]}
FROZEN = {  # module -> constants recorded as the frozen configuration
    "build_pool": ["LEVEL_RANK", "MIN_MEAN_LENGTH", "PLATFORMS"],
    "build_s2": ["CLONE_ANI", "CLONE_AF", "SPECIES_ANI", "SPECIES_AF", "N_SPECIES_NEIGHBOURS", "AUDIT_MIN", "AUDIT_MAX",
                 "DISTANCE_TIERS", "PER_TIER", "SOURCE_RULES"],
    "s3_pair_calls": ["VERSION", "MM", "MIN_CORE", "MAX_CORE", "FLANK", "FLANK_MIN_ID", "FLANK_MIN_COVER", "MERGE_MATCH",
                      "BLOCK_OVERLAP", "SMALL_OTHER_GAP", "EDGE", "COMPETE_COVER", "COMPETE_ID", "BOUNDARY_TOL", "UD",
                      "PROBE_MM", "PROBE_CORE_MAX", "DIFF_ORDER"],
    "s4_annotate": ["VERSION", "SHORT", "CONTEXT", "LISTS", "MMSEQS_SETTINGS", "MEMBER_MIN", "MMSEQS_CHECK", "STARTS", "STOPS"],
    "s3_readcheck": ["VERSION", "FLANK", "ANCHOR", "MIN_DELTA_AS", "MIN_SPAN", "MAX_OTHER_FRACTION", "MIN_INDEL", "LOCUS_TOL",
                     "PRESET", "MM", "MM_GENOME"],
}
BASELINE_POOL = HERE / "comamonadaceae_20261024"


def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(cmd, cwd=None):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, cwd=cwd)
    return r.stdout.strip()


def jsonable(v):
    if isinstance(v, (set, frozenset)):
        return sorted(v)
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    return v


def baseline(a):
    sys.path.insert(0, str(HERE))
    git = dict(branch=run(["git", "rev-parse", "--abbrev-ref", "HEAD"], REPO), head=run(["git", "rev-parse", "HEAD"], REPO),
               origin_head=run(["git", "rev-parse", "@{u}"], REPO), status=run(["git", "status", "--porcelain"], REPO).splitlines())
    tools = {}
    for name, path in TOOLS.items():
        if path.exists():
            out = run([str(path)] + VERSION_ARGS.get(name, ["--version"]))
            tools[name] = dict(path=str(path), version=out.splitlines()[0] if out else "", sha256=sha256(path.resolve()))
        else:
            tools[name] = dict(path=str(path), version="missing")
    params = {}
    for mod, names in FROZEN.items():
        m = importlib.import_module(mod)
        params[mod] = {n: jsonable(getattr(m, n)) for n in names if hasattr(m, n)}
        params[mod]["_source_sha256"] = sha256(HERE / (mod + ".py"))
    params["s4_annotate"]["_translation"] = "table 11 hard-coded (CODE11); transl_table not read yet"
    results = {}
    for p in sorted(BASELINE_POOL.rglob("*")):
        if p.is_file() and p.suffix in (".json", ".tsv", ".gz", ".sha256") and "ena_transport_pilot" not in p.parts:
            results[str(p.relative_to(HERE))] = sha256(p)
    data = {}
    for d in [PROJ / "data/genomes", PROJ / "data/reads"] + sorted((PROJ / "fna_pool").iterdir()):
        if d.is_dir():
            data[str(d)] = dict(entries=sum(1 for _ in d.iterdir()), du=run(["du", "-sh", str(d)]).split("\t")[0])
    out = dict(
        schema="scaleup_stage1_baseline_v1", created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), git=git,
        slurm_queue=run(["squeue", "-u", "kh36969", "-h", "-o", "%i %j %T %M"]).splitlines(),
        running_outside_slurm=a.note or [],
        tools=tools, frozen_parameters=params, baseline_pool=str(BASELINE_POOL.relative_to(REPO)),
        baseline_result_sha256=results, data_directories=data)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "baseline.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(dict(git=git, tools={k: v["version"] for k, v in tools.items()}, n_baseline_files=len(results)), indent=1))


def read_tsv(path):
    import csv
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(path, rows, cols=None):
    import csv
    cols = cols or list(rows[0])
    with open(str(path), "w") as h:
        w = csv.DictWriter(h, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def targets(a):
    from collections import Counter, defaultdict
    rows = [json.loads(l) for l in open(str(a.priorities))]
    wave = sorted((r for r in rows if r["priority"] == "P1_lower_attention_review"), key=lambda r: r["review_rank"])
    src = dict(path=str(a.priorities.resolve().relative_to(REPO)), sha256=sha256(a.priorities))
    out, seen_bs, seen_acc = [], defaultdict(set), defaultdict(set)
    for r in wave:
        fam = r["family_name"].lower()
        d = a.out / fam
        s = json.loads((d / "summary.json").read_text())
        pool = read_tsv(d / "pool.tsv")
        for p in pool:
            seen_bs[p["biosample"]].add(fam)
            for x in (p["genbank_accession"], p["refseq_accession"]):
                if x:
                    seen_acc[x].add(fam)
        n_bs = len({p["biosample"] for p in pool})
        out.append(dict(
            wave=1, review_rank=r["review_rank"], batch_id="family_%d" % r["family_tax_id"], family=r["family_name"],
            taxon_id=r["family_tax_id"], priority=r["priority"], attention_review=r["attention_review"],
            selection_reason="P1_lower_attention_review by archived review_rank; not chosen by candidate function",
            priority_source=src["path"], priority_source_sha256=src["sha256"],
            reads_meta_eligible_wgs_fastq_samples=r["eligible_wgs_fastq_samples"],
            reads_meta_eligible_genomic_wgs_samples=r["eligible_genomic_wgs_samples"], reads_meta_studies=r["distinct_studies"],
            reads_meta_review_flags=";".join(r["review_flags"]),
            ncbi_selection=s["selection"], ncbi_genbank_rows=s["catalog"]["genbank"], ncbi_refseq_rows=s["catalog"]["refseq"],
            pool_assemblies=len(pool), pool_biosamples=n_bs, pool_projects=s["pool"]["n_projects"],
            pool_annotated=s["pool"]["n_annotated_samples"], pool_annotation_coverage=round(s["pool"]["n_annotated_samples"] / len(pool), 4),
            pool_levels=json.dumps(s["pool"]["assembly_level"], sort_keys=True),
            pool_samples_with_usable_long_reads=sum(1 for p in pool if p["linked_read_types"]),
            excluded=json.dumps(s["excluded"], sort_keys=True), min_samples=r["min_samples"],
            meets_min_samples=n_bs >= r["min_samples"],
            batch_status="main" if n_bs >= r["min_samples"] else "below_min_samples_listed_separately",
            ncbi_fetch=s["provenance"]["date"]))
    write_tsv(a.out / "target_batches.tsv", out)
    overlap = [dict(kind="biosample", id=k, batches=";".join(sorted(v))) for k, v in seen_bs.items() if len(v) > 1] + \
              [dict(kind="accession", id=k, batches=";".join(sorted(v))) for k, v in seen_acc.items() if len(v) > 1]
    write_tsv(a.out / "cross_batch_overlap.tsv", overlap or [dict(kind="", id="", batches="")], ["kind", "id", "batches"])
    old = {p["biosample"]: p for p in read_tsv(BASELINE_POOL / "pool.tsv")}
    new = {p["biosample"]: p for p in read_tsv(a.out / "comamonadaceae" / "pool.tsv")}
    excl = {e["biosample"]: e for e in read_tsv(a.out / "comamonadaceae" / "excluded.tsv.gz")}
    newasm = {r["biosample"]: r for r in read_tsv(a.out / "comamonadaceae" / "assemblies.tsv.gz")}
    delta = []
    for bs in sorted(set(old) | set(new)):
        o, n = old.get(bs), new.get(bs)
        if o and n:
            st = "same_assembly" if o["download_accession"] == n["download_accession"] else "different_representative"
        elif o:
            st = "dropped:" + (excl[bs]["reason"] if bs in excl else "not_in_new_ncbi_listing" if bs not in newasm else "other")
        else:
            st = "added"
        delta.append(dict(biosample=bs, status=st, old_role=o["role"] if o else "", old_accession=o["download_accession"] if o else "",
                          new_accession=n["download_accession"] if n else "", level=(n or o)["assembly_level"],
                          linked_read_types=(n or o)["linked_read_types"]))
    write_tsv(a.out / "baseline_pool_delta.tsv", delta)
    print(json.dumps(dict(batches=[(x["family"], x["pool_biosamples"], x["batch_status"]) for x in out],
                          total_pool=sum(x["pool_assemblies"] for x in out), cross_batch_overlap=len(overlap),
                          baseline_delta=dict(Counter(x["status"].split(":")[0] if x["status"].startswith("added") else x["status"] for x in delta))),
                     indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["baseline", "targets"])
    p.add_argument("--out", type=Path, required=True, help="repo output dir of this scale-up version")
    p.add_argument("--note", nargs="*", help="baseline: work running outside SLURM (e.g. tmux sessions on lrc-xfer)")
    p.add_argument("--priorities", type=Path, help="targets: archived batch_priorities.jsonl")
    a = p.parse_args()
    dict(baseline=baseline, targets=targets)[a.mode](a)


if __name__ == "__main__":
    main()
