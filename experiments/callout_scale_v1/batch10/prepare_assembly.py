"""Download GCA and GCF assemblies for one sample and choose the callout input (protocol section 2).

Use the RefSeq FASTA and its annotation only when every record's sequence equals the
GenBank assembly; otherwise use the GenBank FASTA with its own annotation if it has
CDS features, or no annotation (tnp status then unavailable). Never mix a GCF
annotation onto different GCA sequences. Writes assembly_check.json and manifest.json.
"""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from Bio import SeqIO

DATASETS = "/global/home/users/kh36969/.conda/envs/claude-env/bin/datasets"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(acc, out):
    target = out / acc
    if not target.exists():
        subprocess.run([DATASETS, "download", "genome", "accession", acc, "--include", "genome,gbff",
                        "--filename", str(out / (acc + ".zip")), "--no-progressbar"], check=True)
        subprocess.run(["unzip", "-q", "-o", str(out / (acc + ".zip")), "-d", str(target)], check=True)
    data = target / "ncbi_dataset/data" / acc
    fna = sorted(data.glob("*_genomic.fna"))
    gbff = data / "genomic.gbff"
    if len(fna) != 1:
        raise SystemExit("Expected one genomic FASTA for " + acc)
    return fna[0], gbff if gbff.exists() else None


def records(fna):
    return {r.id: (len(r.seq), hashlib.sha256(str(r.seq).upper().encode()).hexdigest()) for r in SeqIO.parse(fna, "fasta")}


def cds_count(gbff):
    if gbff is None:
        return 0
    return sum(1 for r in SeqIO.parse(gbff, "genbank") for f in r.features if f.type == "CDS")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--gca", required=True)
    p.add_argument("--gcf", required=True, help='RefSeq accession, or "none" when there is no pair')
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if (a.out / "manifest.json").exists():
        raise SystemExit("Refusing to overwrite existing " + str(a.out / "manifest.json"))
    gca_fna, gca_gbff = fetch(a.gca, a.out)
    if a.gcf == "none":  # no RefSeq pair: GenBank sequence and its own annotation only
        gcf_fna, gcf_gbff, gcf_rec, identical = None, None, {}, False
    else:
        gcf_fna, gcf_gbff = fetch(a.gcf, a.out)
        gcf_rec = records(gcf_fna)
    gca_rec = records(gca_fna)
    if gcf_fna:
        identical = sorted(v[1] for v in gca_rec.values()) == sorted(v[1] for v in gcf_rec.values())
    if identical and cds_count(gcf_gbff):
        fasta, gbff, choice = gcf_fna, gcf_gbff, "gcf_sequence_and_annotation"
    elif cds_count(gca_gbff):
        fasta, gbff, choice = gca_fna, gca_gbff, "gca_sequence_and_own_annotation"
    else:
        fasta, gbff, choice = gca_fna, None, "gca_sequence_no_annotation_tnp_unavailable"
    check = dict(sample=a.sample, gca=a.gca, gcf=a.gcf, identical_sequence_sets=identical, choice=choice,
                 gca_records=gca_rec, gcf_records=gcf_rec, cds_features=cds_count(gbff),
                 fasta=str(fasta), fasta_sha256=sha256(fasta), gbff=str(gbff) if gbff else None,
                 gbff_sha256=sha256(gbff) if gbff else None)
    (a.out / "assembly_check.json").write_text(json.dumps(check, indent=2) + "\n")
    sample = dict(id=a.sample, fasta=str(fasta))
    if gbff:
        sample.update(annotations=[str(gbff)], annotation_format="genbank")
    (a.out / "manifest.json").write_text(json.dumps(dict(schema_version=1, samples=[sample]), indent=2) + "\n")
    print(json.dumps({k: v for k, v in check.items() if not k.endswith("_records")}, indent=2))


if __name__ == "__main__":
    main()
