"""Census (metadata only) of ENA Oxford Nanopore genomic runs across Bacteria and their exact-BioSample assemblies.

Not a panel: no reads, no assembly FASTA, nothing frozen. Counts how many runs could enter the
assembly + same-BioSample ONT route, by the same first filters as batch50 where metadata allows:
WGS, non-PCR library, exact BioSample with a latest GenBank assembly, Complete Genome or
<= 10 scaffolds (largest record is NOT checked here: no FASTA), and a rough host/source screen.

Inputs: ENA portal read_run search over tax_tree(2); NCBI genomes/genbank/bacteria/assembly_summary.txt
(streamed; only rows whose BioSample appears in the ENA set are kept). Raw files stay in --work.
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import time
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
QUERY = 'tax_tree(2) AND instrument_platform="OXFORD_NANOPORE" AND library_source="GENOMIC"'
FIELDS = ("run_accession,study_accession,sample_accession,tax_id,scientific_name,instrument_model,library_strategy,"
          "library_selection,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp,host,host_scientific_name,"
          "isolation_source,environmental_sample,first_public")
ASM = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
CLINICAL = re.compile(r"homo sapiens|human|patient|sputum|blood|urine|stool|feces|faeces|clinical|hospital|wound|"
                      r"swab|abscess|cystic fibrosis|bovine|cattle|chicken|pig\b|swine|mouse|poultry", re.I)


def fetch(work):
    work.mkdir(parents=True, exist_ok=True)
    url = ENA + "?" + urllib.parse.urlencode(dict(result="read_run", query=QUERY, fields=FIELDS, format="tsv", limit=0))
    subprocess.run(["curl", "-sS", "--fail", "--retry", "3", "-o", str(work / "ena_ont_genomic_runs.tsv"), url], check=True)
    runs = list(csv.DictReader(open(work / "ena_ont_genomic_runs.tsv"), delimiter="\t"))
    wanted = {r["sample_accession"] for r in runs if r["sample_accession"]}
    kept, header, scanned = [], None, 0
    proc = subprocess.Popen(["curl", "-sS", "--fail", "--retry", "3", ASM], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = line[1:].rstrip("\n").split("\t")
            continue
        scanned += 1
        f = line.rstrip("\n").split("\t")
        if len(f) > 2 and f[2] in wanted:
            kept.append(dict(zip(header, f)))
    if proc.wait():
        raise SystemExit("assembly_summary download failed")
    with gzip.open(work / "assembly_summary_matched.tsv.gz", "wt") as h:
        w = csv.DictWriter(h, header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept)
    prov = dict(date=time.strftime("%Y-%m-%d"), ena_url=url, ena_rows=len(runs),
                ena_sha256=hashlib.sha256((work / "ena_ont_genomic_runs.tsv").read_bytes()).hexdigest(),
                assembly_summary_url=ASM, assembly_rows_scanned=scanned, assembly_rows_matched=len(kept))
    (work / "fetch_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    return prov


def census(work, prov):
    runs = list(csv.DictReader(open(work / "ena_ont_genomic_runs.tsv"), delimiter="\t"))
    asm = defaultdict(list)
    with gzip.open(work / "assembly_summary_matched.tsv.gz", "rt") as h:
        for a in csv.DictReader(h, delimiter="\t"):
            if a["version_status"] == "latest":
                asm[a["biosample"]].append(a)
    stage = Counter()
    genus = defaultdict(Counter)
    for r in runs:
        g = (r["scientific_name"].split() or ["unknown"])[0]
        gc = genus[g]
        stage["genomic_runs"] += 1
        gc["genomic_runs"] += 1
        if r["library_strategy"] != "WGS":
            continue
        stage["wgs"] += 1
        if re.search("PCR", r["library_selection"], re.I):
            stage["wgs_pcr"] += 1
            continue
        stage["wgs_non_pcr"] += 1
        gc["wgs_non_pcr"] += 1
        if not (r["fastq_md5"] and r["fastq_bytes"]):
            continue
        stage["with_ena_fastq_md5"] += 1
        matches = asm.get(r["sample_accession"], [])
        if not matches:
            continue
        stage["exact_biosample_assembly"] += 1
        gc["exact_biosample_assembly"] += 1
        levels = {a["assembly_level"] for a in matches}
        for lv in levels:
            stage["level:" + lv] += 1
        contiguous = any(a["assembly_level"] == "Complete Genome" or
                         (a["scaffold_count"].isdigit() and int(a["scaffold_count"]) <= 10) for a in matches)
        if not contiguous:
            continue
        stage["complete_or_le10_records"] += 1
        gc["complete_or_le10_records"] += 1
        if any(a["assembly_level"] == "Complete Genome" for a in matches):
            stage["complete_genome"] += 1
        src = " ".join([r["host"], r["host_scientific_name"], r["isolation_source"]])
        if CLINICAL.search(src):
            stage["host_or_source_clinical_animal_keyword"] += 1
            continue
        stage["after_rough_source_screen"] += 1
        gc["after_rough_source_screen"] += 1
    biosamples = Counter()
    for r in runs:
        if r["library_strategy"] == "WGS" and not re.search("PCR", r["library_selection"], re.I) and asm.get(r["sample_accession"]):
            biosamples["exact_biosample_non_pcr_wgs"] += 1
    summary = dict(schema="bacteria_ont_census_v1", scope="metadata only; not a panel; largest record not checked",
                   provenance=prov, funnel=dict(stage),
                   distinct_biosamples_exact_match_wgs_non_pcr=len({r["sample_accession"] for r in runs
                       if r["library_strategy"] == "WGS" and not re.search("PCR", r["library_selection"], re.I)
                       and asm.get(r["sample_accession"])}),
                   genera_with_exact_match=sum(1 for g in genus.values() if g["exact_biosample_assembly"]),
                   genera_after_rough_screen=sum(1 for g in genus.values() if g["after_rough_source_screen"]))
    (HERE / "census_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    cols = ["genus", "genomic_runs", "wgs_non_pcr", "exact_biosample_assembly", "complete_or_le10_records",
            "after_rough_source_screen"]
    with open(HERE / "genus_counts.tsv", "w") as h:
        h.write("\t".join(cols) + "\n")
        for g, c in sorted(genus.items(), key=lambda x: (-x[1]["after_rough_source_screen"], -x[1]["genomic_runs"], x[0])):
            h.write("\t".join([g] + [str(c[k]) for k in cols[1:]]) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout/bacteria_ont_census"))
    p.add_argument("--reuse", action="store_true", help="skip downloads; use files already in --work")
    a = p.parse_args()
    prov = json.loads((a.work / "fetch_provenance.json").read_text()) if a.reuse else fetch(a.work)
    census(a.work, prov)


if __name__ == "__main__":
    main()
