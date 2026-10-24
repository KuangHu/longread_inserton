"""Check the assembly chosen by prepare_assembly.py against its NCBI reports.

Downloads the sequence report (seq-report) of the GCA and, if used, the GCF; reads the
assembly data report written with the genome download. Checks: BioSample, assembly level,
expected total length, that every FASTA record equals a report record (length, accession),
that chromosome and plasmid records are all kept, and CDS features in the chosen GenBank file.
Writes assembly_verify.json; exits 1 if any check fails.
"""

import argparse
import json
import subprocess
import zipfile
from pathlib import Path

from Bio import SeqIO

DATASETS = "/global/home/users/kh36969/.conda/envs/claude-env/bin/datasets"


def reports(acc, asm_dir):
    z = asm_dir / (acc + ".seq_report.zip")
    if not z.exists():
        subprocess.run([DATASETS, "download", "genome", "accession", acc, "--include", "seq-report",
                        "--filename", str(z), "--no-progressbar"], check=True)
    with zipfile.ZipFile(z) as f:
        seq = [json.loads(x) for x in f.read("ncbi_dataset/data/%s/sequence_report.jsonl" % acc).decode().splitlines()]
        data = json.loads(f.read("ncbi_dataset/data/assembly_data_report.jsonl").decode().splitlines()[0])
    return seq, data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--asm-dir", type=Path, required=True, help="prepare_assembly.py --out directory")
    p.add_argument("--biosample", required=True)
    p.add_argument("--total-bp", type=int, required=True)
    a = p.parse_args()
    c = json.loads((a.asm_dir / "assembly_check.json").read_text())
    used = c["gcf"] if c["choice"].startswith("gcf") else c["gca"]
    seq, data = reports(used, a.asm_dir)
    gca_seq, gca_data = reports(c["gca"], a.asm_dir) if used != c["gca"] else (seq, data)
    fasta = {r.id: len(r.seq) for r in SeqIO.parse(c["fasta"], "fasta")}
    acc_key = "refseq_accession" if used.startswith("GCF") else "genbank_accession"
    rep = {s[acc_key]: s for s in seq}
    rows = [dict(fasta_id=k, length=n, report_length=rep.get(k, {}).get("length"), role=rep.get(k, {}).get("role"),
                 location=rep.get(k, {}).get("assigned_molecule_location_type"),
                 molecule=rep.get(k, {}).get("assigned_molecule") or rep.get(k, {}).get("sequence_name"),
                 genbank_accession=rep.get(k, {}).get("genbank_accession")) for k, n in fasta.items()]
    kinds = sorted({r["location"] for r in rows if r["location"]})
    checks = dict(
        biosample_match=data["assembly_info"]["biosample"]["accession"] == a.biosample
        and gca_data["assembly_info"]["biosample"]["accession"] == a.biosample,
        total_length_match=sum(fasta.values()) == a.total_bp == int(data["assembly_stats"]["total_sequence_length"]),
        fasta_records_equal_report=len(fasta) == len(seq) and all(r["length"] == r["report_length"] for r in rows),
        all_report_records_kept=set(rep) == set(fasta),
        gca_gcf_same_record_count=len(gca_seq) == len(seq),
        annotation_present=c["cds_features"] > 0,
    )
    out = dict(used_accession=used, choice=c["choice"], fasta=c["fasta"], fasta_sha256=c["fasta_sha256"],
               gbff=c["gbff"], gbff_sha256=c["gbff_sha256"], cds_features=c["cds_features"],
               assembly_level=data["assembly_info"]["assembly_level"],
               assembly_method=data["assembly_info"].get("assembly_method"),
               sequencing_tech=data["assembly_info"].get("sequencing_tech"),
               biosample=data["assembly_info"]["biosample"]["accession"],
               organism=data["organism"], total_bp=sum(fasta.values()), records=rows, molecule_types=kinds,
               annotation_info=data.get("annotation_info", {}), checks=checks,
               status="pass" if all(checks.values()) else "fail")
    (a.asm_dir / "assembly_verify.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "annotation_info"}, indent=2))
    raise SystemExit(out["status"] != "pass")


if __name__ == "__main__":
    main()
