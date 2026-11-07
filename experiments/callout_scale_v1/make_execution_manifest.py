"""Resolve test_panel.tsv into an execution manifest of already-verified inputs.

test_panel.tsv is the acceptance list; it is not the read_callout samples.tsv
format. Each row is mapped to the verified reads, assembly FASTA/annotation and
existing v1 run on the cluster. Assembly and annotation SHA256 are recomputed
and must equal the fetch-time records; reads are checked for size and verify
status here and re-hashed inside the job. Nothing is downloaded.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Facts recorded in the panel audits (ENA library_selection, NCBI assembly method).
LIBRARY = {
    "PaAR2038": ("RANDOM PCR", "Illumina NovaSeq + Nanopore PromethION; Unicycler/Racon/Pilon"),
    "SMT4a": ("PCR", "Illumina HiSeq + ONT MinION; SPAdes/Unicycler"),
    "LjRoot223": ("RANDOM", "Illumina NovaSeq 6000 + MinION; Unicycler"),
    "R-93128": ("RANDOM", "ONT MinION only; dragonflye"),
    "R-91497": ("RANDOM", "ONT MinION only; dragonflye"),
}
V1_RUNS = {"PaAR2038": "runs/PaAR2038"}  # others are in runs_paraburkholderia_round1/


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proj", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout"))
    p.add_argument("--out", type=Path, default=HERE / "execution_manifest.json")
    a = p.parse_args()
    rows = list(csv.DictReader(open(HERE / "test_panel.tsv"), delimiter="\t"))
    entries = []
    for r in rows:
        sid = r["sample_id"]
        asm_dir = a.proj / "data/assemblies" / sid
        check = json.loads((asm_dir / "assembly_check.json").read_text())
        assert check["gca"] == r["assembly"], sid + ": panel GCA differs from fetched GCA"
        assert check["identical_sequence_sets"], sid + ": GCA/GCF sequence sets differ"
        fasta, gbff = Path(check["fasta"]), Path(check["gbff"])
        assert sha256(fasta) == check["fasta_sha256"], sid + ": FASTA hash changed since fetch"
        assert sha256(gbff) == check["gbff_sha256"], sid + ": annotation hash changed since fetch"
        verify = json.loads((a.proj / "data/reads" / r["run"] / "verify.json").read_text())
        assert verify["status"] == "verified" and verify["run"]["run"] == r["run"], sid + ": reads not verified"
        assert verify["run"]["biosample"] == r["biosample"], sid + ": reads BioSample differs from panel"
        reads = Path(verify["local_file"])
        assert reads.stat().st_size == verify["bytes"], sid + ": reads file size changed"
        v1 = a.proj / V1_RUNS.get(sid, "runs_paraburkholderia_round1/" + sid)
        selection, assembly = LIBRARY[sid]
        pcr = "PCR" in selection
        entries.append(dict(
            order=int(r["order"]), sample_id=sid, role=r["role"], run=r["run"], biosample=r["biosample"],
            gca=check["gca"], gcf_used_for_sequence_and_annotation=check["gcf"],
            fasta=str(fasta), fasta_sha256=check["fasta_sha256"], gbff=str(gbff), gbff_sha256=check["gbff_sha256"],
            reads=str(reads), reads_bytes=verify["bytes"], reads_md5=verify["md5"], reads_sha256=verify["sha256"],
            library_selection=selection, assembly_method=assembly, platform=r["platform"], preset="map-ont",
            v1_sample_dir=str(v1), v1_callout=str(v1 / "callout"), v1_evidence=str(v1 / "evidence"),
            independence="pcr_library_independence_unverified" if pcr else "independence_unverified",
            assembly_consistency="assembly_built_from_these_reads" if "only" in assembly else "hybrid_assembly_used_these_reads",
        ))
    out = dict(schema_version=1, source="experiments/callout_scale_v1/test_panel.tsv",
               note="Same verified files as the v1 runs; no re-download. GCF sequences equal GCA (per-record SHA256).",
               samples=entries)
    a.out.write_text(json.dumps(out, indent=2) + "\n")
    for e in entries:
        print(e["order"], e["sample_id"], e["run"], e["gcf_used_for_sequence_and_annotation"], e["library_selection"], e["independence"])


if __name__ == "__main__":
    main()
