"""Compare a package run with the frozen three-isolate pilot, without rerunning discovery."""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    expected = json.loads(Path(__file__).with_name("pilot_expected.json").read_text())
    observed = json.loads((args.run / "results/member_boundaries.json").read_text())
    keys = expected["coordinate_columns"]
    actual = [[r[k] for k in keys] for r in observed]
    if actual != expected["member_coordinates"]:
        raise SystemExit("FAIL: pilot coordinates/statuses changed")
    for rel, sha in expected["sequence_sha256"].items():
        if hashlib.sha256((args.run / rel).read_bytes()).hexdigest() != sha:
            raise SystemExit("FAIL: sequence export changed: " + rel)
    print(
        "PASS: all 300 native/inferred/exported coordinate rows and eight FASTA files match the frozen pilot"
    )


if __name__ == "__main__":
    main()
