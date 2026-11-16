"""Member-by-member difference between two audited_filled_site_v1 outputs of one sample.

Lists every member whose audited class or unit counts changed; an empty table is
still written so that "no change" is recorded rather than assumed.
"""

import argparse
import csv
import json
from pathlib import Path

COLS = ["audited_class", "units_full_support", "units_x_deleted", "units_competing", "units_cigar_fail"]


def load(d):
    return {r["member"]: r for r in csv.DictReader(open(Path(d) / "member_audit.tsv"), delimiter="\t")}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old", type=Path, required=True)
    p.add_argument("--new", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    old, new = load(a.old), load(a.new)
    if set(old) != set(new):
        raise SystemExit("Member sets differ between old and new audits")
    rows = []
    for m in sorted(new):
        if any(old[m][c] != new[m][c] for c in COLS):
            rows.append(dict(member=m, length=new[m]["length"], **{"old_" + c: old[m][c] for c in COLS},
                             **{"new_" + c: new[m][c] for c in COLS},
                             class_changed=old[m]["audited_class"] != new[m]["audited_class"]))
    cols = ["member", "length", "class_changed"] + ["old_" + c for c in COLS] + ["new_" + c for c in COLS]
    with open(a.out, "w") as h:
        h.write("\t".join(cols) + "\n")
        for r in rows:
            h.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(json.dumps(dict(members=len(new), changed_rows=len(rows), class_changed=sum(r["class_changed"] for r in rows))))


if __name__ == "__main__":
    main()
