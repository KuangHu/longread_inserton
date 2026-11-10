"""Turn the frozen selection_manifest.tsv into the batch50 job table and register its runs.

Does not change the selection. Per row: batch10 columns (order..preset) plus panel rank, ENA
bytes/MD5/bases, the FASTA the job must end up using (RefSeq when its record set equals GenBank,
as prepare_assembly.py decides) with its expected SHA256 from the metadata audit download, and the
read policy: ENA bases / assembly bp above DOWNSAMPLE_ABOVE -> fixed-seed downsample to TARGET_DEPTH
for mapping and the read check only. Appends role callout_scale_v1_batch50 entries to
scripts/long_read_runs.json (refuses to change an existing entry).
"""

import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY = HERE.parents[2] / "scripts" / "long_read_runs.json"
DOWNSAMPLE_ABOVE, TARGET_DEPTH, SEED = 200.0, 80, 20261024
COLS = ["order", "sample_id", "run", "biosample", "gca", "gcf", "preset", "panel_rank", "organism", "continuity",
        "fastq_bytes", "fastq_md5", "ena_read_count", "ena_base_count", "assembly_bp", "ena_depth",
        "selected_fasta_source", "selected_fasta_accession", "expected_fasta_sha256", "expected_records",
        "read_policy", "downsample_depth", "downsample_seed", "fastq_url"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genomes", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout/batch50_meta/genomes"))
    a = p.parse_args()
    sel = list(csv.DictReader(open(HERE / "selection_manifest.tsv"), delimiter="\t"))
    out, reg_new = [], []
    for r in sel:
        identical = r["gca_gcf_check"] == "done_identical"
        if not identical and r["gca_gcf_check"] != "no_gcf_use_gca":
            raise SystemExit("%s: GCA/GCF check not done (%s)" % (r["sample_id"], r["gca_gcf_check"]))
        acc = r["gcf"] if identical else r["gca"]
        fa = json.loads((a.genomes / (acc + ".json")).read_text())
        bp = int(r["total_bp"])
        if fa["total_bp"] != bp:
            raise SystemExit("%s: selected FASTA length differs from GCA" % r["sample_id"])
        depth = int(r["base_count"]) / bp
        ds = depth > DOWNSAMPLE_ABOVE
        out.append(dict(
            order=r["order"], sample_id=r["sample_id"], run=r["run"], biosample=r["biosample"], gca=r["gca"],
            gcf=r["gcf"] or "none", preset="map-ont", panel_rank=r["rank"], organism=r["organism"],
            continuity=r["continuity"], fastq_bytes=r["fastq_bytes"], fastq_md5=r["fastq_md5"],
            ena_read_count=r["read_count"], ena_base_count=r["base_count"], assembly_bp=bp, ena_depth="%.1f" % depth,
            selected_fasta_source="gcf_sequence_and_annotation" if identical else "gca_sequence",
            selected_fasta_accession=acc, expected_fasta_sha256=fa["fasta_sha256"], expected_records=fa["records"],
            read_policy="downsample" if ds else "full", downsample_depth=TARGET_DEPTH if ds else "",
            downsample_seed=SEED if ds else "", fastq_url="https://" + r["fastq_ftp"]))
        reg_new.append(dict(
            sample_id=r["sample_id"], run=r["run"], biosample=r["biosample"], assembly=r["gca"],
            platform="OXFORD_NANOPORE", instrument_model=r["instrument_model"], role="callout_scale_v1_batch50",
            url="https://" + r["fastq_ftp"], archive_read_count=int(r["read_count"]),
            archive_base_count=int(r["base_count"]), fastq_bytes=int(r["fastq_bytes"]), fastq_md5=r["fastq_md5"]))
    with open(HERE / "samples.tsv", "w") as h:
        w = csv.DictWriter(h, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    reg = json.loads(REGISTRY.read_text())
    have = {e["run"]: e for e in reg["runs"]}
    added = 0
    for e in reg_new:
        if e["run"] in have:
            if have[e["run"]] != e:
                raise SystemExit("Registry already has a different entry for " + e["run"])
            continue
        if any(x["sample_id"] == e["sample_id"] for x in reg["runs"]):
            raise SystemExit("Registry sample_id clash: " + e["sample_id"])
        reg["runs"].append(e)
        added += 1
    if added:
        reg["source"] += "; batch50 from experiments/callout_scale_v1/batch50/selection_manifest.tsv (ENA 2026-10-24)"
        REGISTRY.write_text(json.dumps(reg, indent=2) + "\n")
    print(len(out), "samples;", sum(o["read_policy"] == "downsample" for o in out), "downsample;", added, "registry entries added")


if __name__ == "__main__":
    main()
