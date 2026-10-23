"""Check an audited_filled_site_v1 output directory against its call-out run.

Independent of the classifier: recomputes each member's support units from
read_evidence.tsv, re-derives the audited class with the same precedence, and
checks member identity/coordinates against member_boundaries.json. It does not
re-map reads and cannot certify independence of molecules.
"""

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

ORDER = ["full_support", "x_deleted", "competing", "full_span_cigar_fail", "left_only", "right_only", "internal"]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def expected_class(row, n, min_units, contig_len):
    s, e = int(row["start0"]), int(row["end0"])
    if s - 100 < 0 or e + 100 > contig_len:
        return "edge_limited"
    if n["full_support"] >= min_units:
        return "audited_full_span_with_conflict" if n["x_deleted"] >= min_units else "audited_full_span"
    if n["full_support"] == 1:
        return "single_unit_full_span"
    for status, cls in (("competing", "ambiguous_competing"), ("full_span_cigar_fail", "full_span_cigar_fail"),
                        ("x_deleted", "conflict_x_absent")):
        if n[status]:
            return cls
    if n["left_only"] or n["right_only"]:
        return "one_side_only"
    return "internal_only" if n["internal"] else "no_coverage"


def check(audit_dir, members_json, reference):
    audit_dir = Path(audit_dir)
    summary = json.loads((audit_dir / "summary.json").read_text())
    require(summary["schema"] == "audited_filled_site_v1", "Not an audited_filled_site_v1 directory")
    opener = gzip.open if str(members_json).endswith(".gz") else open
    with opener(members_json, "rt") as h:
        boundaries = {r["member"]: r for r in json.load(h) if r.get("candidate_start0") is not None}
    lengths, name = {}, None
    for line in open(reference):
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        else:
            lengths[name] += len(line.strip())
    audit = list(csv.DictReader(open(audit_dir / "member_audit.tsv"), delimiter="\t"))
    by_member = {r["member"]: r for r in audit}
    require(len(by_member) == len(audit), "Duplicate members in member_audit.tsv")
    require(set(by_member) == set(boundaries), "member_audit.tsv members differ from member_boundaries.json")
    verdicts = defaultdict(dict)
    for r in csv.DictReader(open(audit_dir / "read_evidence.tsv"), delimiter="\t"):
        require(r["member"] in by_member, "Evidence for unknown member " + r["member"])
        require(r["status"] in ORDER, "Unknown read status " + r["status"])
        if r["duplicate_of_unit"]:
            continue
        old = verdicts[r["member"]].get(r["read"])
        if old is None or ORDER.index(r["status"]) < ORDER.index(old):
            verdicts[r["member"]][r["read"]] = r["status"]
    min_units = summary["thresholds"]["MIN_UNITS"]
    for member, row in by_member.items():
        b = boundaries[member]
        require((row["contig"], int(row["start0"]), int(row["end0"]), row["boundary_source"]) ==
                (b["accession"], b["candidate_start0"], b["candidate_end0"], b["boundary_source"]),
                member + ": coordinates or boundary differ from call-out")
        n = Counter(verdicts.get(member, {}).values())
        for status, col in (("full_support", "units_full_support"), ("x_deleted", "units_x_deleted"),
                            ("competing", "units_competing"), ("full_span_cigar_fail", "units_cigar_fail")):
            require(n[status] == int(row[col]), "%s: %s recomputed %d != %s" % (member, col, n[status], row[col]))
        cls = expected_class(row, n, min_units, lengths[row["contig"]])
        require(cls == row["audited_class"], "%s: class %s != recomputed %s" % (member, row["audited_class"], cls))
    counts = Counter(r["audited_class"] for r in audit)
    require(all(summary["audited_class_counts"].get(c, 0) == v for c, v in counts.items()), "Summary class counts mismatch")
    require(summary["members"] == len(audit), "Summary member count mismatch")
    return dict(schema="audited_evidence_check_v1", sample=summary["sample"], status="pass", members=len(audit),
                audited_class_counts=dict(counts), scope="Recomputed units/classes from read_evidence.tsv; no re-mapping")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-dir", type=Path, required=True)
    p.add_argument("--members", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    try:
        require(not a.out.exists(), "Output already exists; choose a new report path")
        report = check(a.audit_dir, a.members, a.reference)
        with a.out.open("x") as h:
            json.dump(report, h, indent=2)
            h.write("\n")
        print(json.dumps(report))
    except (ValueError, OSError, KeyError) as error:
        p.exit(1, "Error: %s\n" % error)


if __name__ == "__main__":
    main()
