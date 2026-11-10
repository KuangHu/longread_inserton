"""Cross-sample table: what each data type contributed to the insertion callout.

Reads per-sample outputs under PROJ/runs/<sample>/ and writes comparison.tsv
and comparison.json next to them. Samples differ in organism as well as
platform, so differences are not a platform benchmark.
"""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text()) if path.exists() else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proj", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout"))
    args = parser.parse_args()
    samples = list(csv.DictReader(open(HERE / "samples.tsv"), delimiter="\t"))
    rows = []
    for s in samples:
        out = args.proj / "runs" / s["sample_id"]
        callout = load(out / "callout/results/summary.json") or {}
        ev = load(out / "evidence/summary.json") or {}
        unit = load(out / "read_unit.json") or {}
        rename = load(out / "reads/rename_stats.json") or {}
        verify = load(args.proj / "data/reads" / s["run"] / "verify.json") or {}
        support = out / "evidence/member_read_support.tsv"
        by_method = Counter()
        if support.exists():
            for r in csv.DictReader(open(support), delimiter="\t"):
                by_method[(r["method"], r["class"])] += 1
        classes = ev.get("member_classes", {})
        rows.append(dict(
            sample=s["sample_id"], platform=s["platform"], preset=s["preset"], run=s["run"],
            reads_verified=verify.get("status"), molecule_grouping=verify.get("summary", {}).get("molecule_grouping"),
            records_with_adapter=unit.get("fraction_with_adapter"),
            adapter_split_records=rename.get("records_split"),
            molecule_depth=ev.get("molecule_depth"),
            native_families=callout.get("native_families"), members=ev.get("members"),
            boundary_sources=callout.get("boundary_sources"),
            **{"class_" + k: v for k, v in classes.items()},
            nonreference_family_sites=ev.get("nonreference_family_sites_called"),
            nonreference_family_two_sided=ev.get("nonreference_family_sites_two_sided"),
            unassigned_insertion_sites=ev.get("unassigned_insertion_sites_called"),
            palindromic_read_pieces=ev.get("palindromic_read_pieces"),
            classes_by_method={"%s:%s" % k: v for k, v in sorted(by_method.items())},
        ))
    (args.proj / "comparison.json").write_text(json.dumps(rows, indent=2) + "\n")
    flat = [k for k in rows[0] if not isinstance(rows[0][k], dict)] if rows else []
    with open(args.proj / "comparison.tsv", "w") as handle:
        handle.write("field\t" + "\t".join(r["sample"] for r in rows) + "\n")
        for k in flat[1:]:
            handle.write(k + "\t" + "\t".join(str(r.get(k)) for r in rows) + "\n")
    print((args.proj / "comparison.tsv").read_text())


if __name__ == "__main__":
    main()
