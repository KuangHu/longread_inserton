"""Independent checks of the frozen 50-sample metadata panel; re-reads the raw ENA/NCBI/FASTA caches.

Fails (exit 1) on: duplicate run/BioSample/assembly in the selection, any already-run sample,
a primary row that is PCR, not WGS/GENOMIC, not exact-BioSample matched, below continuity,
clinical/animal, B. pseudomallei/mallei or without ENA bytes/MD5; a clone group above the cap;
unstable headers or duplicate row keys in any TSV.
"""

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

from build_panel import ALREADY_RUN, CLONE_CAP, COLS, MAX_RECORDS, MIN_LARGEST, load_tsv

HERE = Path(__file__).resolve().parent
HEADERS = {
    "selection_manifest.tsv": (["order"] + COLS, ("run",)),
    "ranked_candidates.tsv": (COLS + ["selection_pass"], ("run",)),
    "clone_reserve.tsv": (COLS, ("run",)),
    "exclusion_reasons.tsv": (["run", "biosample", "gca", "organism", "assembly_level", "records", "largest_bp",
                               "library_strategy", "library_selection", "source_class", "source_note", "status",
                               "reason"], ("run", "gca")),
    "clone_groups.tsv": (["clone_group", "species_cluster", "ani95_group", "gca", "organism", "sample_id", "run", "status",
                          "group_size", "group_selected", "group_already_run"], ("gca",)),
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout/batch50_meta"))
    p.add_argument("--dir", type=Path, default=HERE)
    a = p.parse_args()
    errors, checks = [], Counter()

    def fail(msg):
        errors.append(msg)

    for name, (cols, key) in HEADERS.items():
        with open(a.dir / name) as h:
            head = h.readline().rstrip("\n").split("\t")
        if head != cols:
            fail("%s: header differs from the fixed column list" % name)
        rows = load_tsv(a.dir / name)
        dup = [k for k, n in Counter(tuple(r[c] for c in key) for r in rows).items() if n > 1]
        if dup:
            fail("%s: duplicate keys %s" % (name, dup[:5]))
        checks["tsv_checked"] += 1

    sel = load_tsv(a.dir / "selection_manifest.tsv")
    for col in ("run", "biosample", "gca", "sample_id"):
        dup = [k for k, n in Counter(r[col] for r in sel).items() if n > 1]
        if dup:
            fail("selection: duplicate %s %s" % (col, dup))

    done = set()
    for ref in ALREADY_RUN.values():
        for r in load_tsv(ref):
            done |= {r["sample_id"], r["run"], r["biosample"], r.get("assembly") or r.get("gca")}

    ena = {r["run_accession"]: r for r in load_tsv(a.work / "raw/ena_ont_runs.tsv")}
    ncbi = {}
    for line in open(a.work / "raw/ncbi_biosample_reports.jsonl"):
        r = json.loads(line)
        ncbi[r["accession"]] = r
    for r in sel:
        tag = r["run"]
        if done & {r["sample_id"], r["run"], r["biosample"], r["gca"]}:
            fail("%s: already-run sample in selection" % tag)
        e, n = ena.get(r["run"]), ncbi.get(r["gca"])
        if not e or not n:
            fail("%s: run or assembly missing from raw caches" % tag)
            continue
        if not (e["library_strategy"] == "WGS" and e["library_source"] == "GENOMIC"
                and e["instrument_platform"] == "OXFORD_NANOPORE"):
            fail("%s: not ONT WGS/GENOMIC" % tag)
        if re.search("PCR", e["library_selection"], re.I):
            fail("%s: PCR library in primary panel" % tag)
        if e["sample_accession"] != n["assembly_info"]["biosample"]["accession"] or e["sample_accession"] != r["biosample"]:
            fail("%s: BioSample differs between run and assembly" % tag)
        if not (e["fastq_md5"] and e["fastq_bytes"] and ";" not in e["fastq_ftp"]):
            fail("%s: ENA bytes/MD5 missing or multi-file" % tag)
        if re.match(r"Burkholderia (pseudo)?mallei\b", n["organism"]["organism_name"]):
            fail("%s: excluded species" % tag)
        f = json.loads((a.work / "genomes" / (r["gca"] + ".json")).read_text())
        complete = n["assembly_info"]["assembly_level"] == "Complete Genome"
        if not (complete or (f["records"] <= MAX_RECORDS and f["largest_bp"] >= MIN_LARGEST)):
            fail("%s: continuity fails" % tag)
        if str(f["records"]) != r["records"] or str(f["largest_bp"]) != r["largest_bp"]:
            fail("%s: FASTA stats differ from manifest" % tag)
        if r["source_class"] == "clinical/animal":
            fail("%s: clinical/animal source" % tag)
        checks["selected_rows_checked"] += 1

    groups = load_tsv(a.dir / "clone_groups.tsv")
    # The cap applies to groups the panel adds to; already-run members count against it.
    per_group = Counter(g["clone_group"] for g in groups if g["status"] in ("selected", "already_run"))
    touched = {g["clone_group"] for g in groups if g["status"] == "selected"}
    over = {c: n for c, n in per_group.items() if c in touched and n > CLONE_CAP}
    if over:
        fail("clone groups above cap: %s" % over)
    # Within a clone group only one sample (selected or already run) may be unmarked; the rest must be
    # explicitly marked second_allowed.
    prior = Counter(g["clone_group"] for g in groups if g["status"] == "already_run")
    firsts = Counter(r["clone_group"] for r in sel if r["clone_member"] == "first")
    for c in {r["clone_group"] for r in sel}:
        if firsts[c] + min(prior[c], 1) != 1:
            fail("clone group %s: %d 'first' selected with %d already run" % (c, firsts[c], prior[c]))

    out = dict(status="pass" if not errors else "fail", errors=errors, checks=dict(checks), selected=len(sel))
    (a.dir / "panel_check.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
