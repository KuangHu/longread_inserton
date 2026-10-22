"""Generate a deterministic repeat fixture; no biological activity is implied."""

import argparse
import json
import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        parser.error("Use a new or empty directory")
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(731)
    dna = lambda n: "".join(rng.choices("ACGT", k=n))
    repeat = dna(1200)
    records = []
    for i in range(4):
        seq = Seq(dna(2500) + repeat + dna(2500))
        record = SeqRecord(
            seq,
            id=f"fixture{i}",
            description="synthetic random repeat; no natural mobile element",
        )
        record.annotations["molecule_type"] = "DNA"
        # A deliberately artificial annotation label exercises post-hoc tnp support.
        location = SimpleLocation(2600, 3500, strand=1)
        record.features = [
            SeqFeature(
                location,
                type="CDS",
                qualifiers={
                    "product": [
                        "synthetic transposase LABEL for software testing only"
                    ],
                    "translation": [str(seq[2600:3500].translate())],
                },
            )
        ]
        records.append(record)
    SeqIO.write(records, args.out / "genome.fasta", "fasta")
    SeqIO.write(records, args.out / "annotation.gb", "genbank")
    (args.out / "manifest.json").write_text(
        json.dumps(
            dict(
                schema_version=1,
                samples=[
                    dict(
                        id="demo",
                        fasta="genome.fasta",
                        annotations=["annotation.gb"],
                        annotation_format="genbank",
                    )
                ],
            ),
            indent=2,
        )
        + "\n"
    )
    print(args.out / "manifest.json")


if __name__ == "__main__":
    main()
