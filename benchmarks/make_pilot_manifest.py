"""Prepare the regression manifest from the original environmental_pilot archive."""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot-root",
        type=Path,
        required=True,
        help="Extracted environmental_pilot directory",
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.pilot_root.resolve()
    source = json.loads(Path(__file__).with_name("source_manifest.json").read_text())
    frozen = json.loads(Path(__file__).with_name("regression_result.json").read_text())
    samples = []
    for sid in ["AltDE", "UT26S", "RHA1"]:
        fasta = root / "inputs" / (sid + ".fasta")
        if (
            hashlib.sha256(fasta.read_bytes()).hexdigest()
            != frozen["input_fasta_sha256"][sid]
        ):
            raise SystemExit("Input FASTA checksum mismatch: " + sid)
        annotations = []
        for row in source:
            if row["organism"] == sid and row["format"] == "embl":
                path = root / "reference" / (row["accession"] + ".embl")
                if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
                    raise SystemExit("Annotation checksum mismatch: " + str(path))
                annotations.append(str(path))
        samples.append(
            dict(
                id=sid,
                fasta=str(fasta),
                annotations=annotations,
                native={
                    method: str(root / "runs" / f"{sid}_{method}_default")
                    for method in ["PILER", "RS"]
                },
            )
        )
    if args.out.exists():
        parser.error("Manifest already exists; choose a new output file")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(dict(schema_version=1, samples=samples), indent=2) + "\n"
    )
    print(args.out)


if __name__ == "__main__":
    main()
