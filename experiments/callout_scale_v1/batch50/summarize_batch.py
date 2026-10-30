"""Archive small per-sample outputs of batch50 jobs and tabulate them (same layout as batch10/results).

For every sample directory under $PROJ/batch50 (or --samples): copy stages, source commit, callout
summary/BED/families, audited summary and member_audit, QC, mapping command/versions, assembly check,
reads verify and downsample record into results/<sample>/; large files are listed with SHA256 in
omitted_large_files.tsv. Writes batch_summary.tsv (one row per sample) and resources.tsv (per stage).
Counts are records, not independent insertions.
"""

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJ = Path("/global/scratch/users/kh36969/longread_callout")
SMALL = {"callout/results/summary.json": "callout/summary.json", "callout/results/candidates.bed": "callout/candidates.bed",
         "callout/results/family_summary.txt": "callout/family_summary.txt",
         "callout/results/all_families_annotated.json": "callout/all_families_annotated.json",
         "evidence_audited/summary.json": "evidence_audited/summary.json",
         "evidence_audited/member_audit.tsv": "evidence_audited/member_audit.tsv",
         "qc/validate.json": "qc/validate.json", "qc/audited_check.json": "qc/audited_check.json",
         "qc/disk_usage.txt": "qc/disk_usage.txt", "mapping/command.txt": "qc/command.txt",
         "mapping/minimap2_version.txt": "qc/minimap2_version.txt", "mapping/cap_reads.txt": "qc/cap_reads.txt",
         "mapping/downsample.json": "qc/downsample.json", "stages.tsv": "stages.tsv",
         "source_commit.txt": "source_commit.txt"}
LARGE = ["callout/results/member_boundaries.json", "evidence_audited/read_evidence.tsv"]
OTHER = ["edge_limited", "ambiguous_competing", "one_side_only"]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", nargs="*")
    p.add_argument("--out", type=Path, default=HERE / "results")
    a = p.parse_args()
    table = {r["sample_id"]: r for r in csv.DictReader(open(HERE / "samples.tsv"), delimiter="\t")}
    names = a.samples or sorted(d.name for d in (PROJ / "batch50").iterdir() if d.is_dir())
    rows, res, omitted = [], [], []
    for s in sorted(names, key=lambda s: int(table[s]["order"])):
        t, d, o = table[s], PROJ / "batch50" / s, a.out / s
        for src, dst in SMALL.items():
            if (d / src).exists():
                (o / dst).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(d / src, o / dst)
        asm = PROJ / "data/assemblies_batch50" / s / "assembly_check.json"
        ver = PROJ / "data/reads" / t["run"] / "verify.json"
        for src, dst in ((asm, "assembly_check.json"), (ver, "reads_verify.json")):
            if src.exists():
                shutil.copy2(src, o / dst)
        for f in LARGE:
            if (d / f).exists():
                omitted.append(dict(sample=s, file_on_cluster="$PROJ/batch50/%s/%s" % (s, f), sha256=sha256(d / f),
                                    bytes=(d / f).stat().st_size))
        stages = list(csv.DictReader(open(d / "stages.tsv"), delimiter="\t")) if (d / "stages.tsv").exists() else []
        res += [dict(sample=s, stage=x["stage"], exit_code=x["exit_code"], wall=x["wall"], max_rss_kb=x["max_rss_kb"])
                for x in stages]
        load = lambda f: json.loads((d / f).read_text()) if (d / f).exists() else {}
        cs, au, va, ac = load("callout/results/summary.json"), load("evidence_audited/summary.json"), \
            load("qc/validate.json"), load("qc/audited_check.json")
        fams = load("callout/results/all_families_annotated.json") or []
        ds = load("mapping/downsample.json")
        cls, lng = au.get("audited_class_counts", {}), au.get("long_member_class_counts", {})
        full = lambda c: c.get("audited_full_span", 0) + c.get("audited_full_span_with_conflict", 0)
        rows.append(dict(
            order=t["order"], sample=s, read_policy=t["read_policy"],
            assembly_choice=json.loads(asm.read_text())["choice"] if asm.exists() else "",
            native_families=cs.get("native_families", ""),
            tnp_supported_families=sum(1 for f in fams if f["tnp_supported"]) if fams else "",
            exported_records=cs.get("exported_records", ""), long_records=sum(lng.values()) if lng else "",
            audited_full=full(cls) if cls else "", long_audited_full=full(lng) if lng else "",
            **{k: cls.get(k, "") for k in OTHER},
            other=sum(v for k, v in cls.items() if k not in OTHER + ["audited_full_span", "audited_full_span_with_conflict"]) if cls else "",
            with_conflict=cls.get("audited_full_span_with_conflict", ""), cap_reads=au.get("cap_reads", ""),
            map_reads_depth=ds.get("output_depth", t["ena_depth"]) if stages else "",
            validate=va.get("status", ""), audited_check=ac.get("status", ""),
            classifier_revision=au.get("classifier_revision", ""),
            all_stages_exit0=bool(stages) and all(x["exit_code"] == "0" for x in stages) and stages[-1]["stage"] == "audited_check"))
    a.out.mkdir(parents=True, exist_ok=True)
    for name, data in (("batch_summary.tsv", rows), ("resources.tsv", res), ("omitted_large_files.tsv", omitted)):
        if data:
            with open(a.out / name, "w") as h:
                w = csv.DictWriter(h, list(data[0]), delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(data)
    for r in rows:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main()
