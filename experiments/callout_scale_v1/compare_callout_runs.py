"""Member-by-member comparison of two multisite-insertion runs on the same frozen input.

Used for the AR20-38 regression: native families, every member's native,
inferred and candidate coordinates, boundary source and flank state, and the
exported X/triplet/flank sequences must be identical. Differences are listed,
not hidden; any difference makes the exit code 1.
"""

import argparse
import json
import sys
from pathlib import Path

FIELDS = ["family", "method", "accession", "strand", "native_start0", "native_end0", "inferred_start0", "inferred_end0",
          "candidate_start0", "candidate_end0", "boundary_source", "flanks_complete", "upstream_length", "downstream_length"]
FASTAS = ["candidate_insertion.fasta", "upstream_500.fasta", "downstream_500.fasta", "triplet.fasta"]


def fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = [line.strip()]
        else:
            seqs[name].append(line.strip())
    return {k: "\n".join(v) for k, v in seqs.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old", type=Path, required=True)
    p.add_argument("--new", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    load = lambda run, f: json.loads((run / "results" / f).read_text())
    diffs = []
    old_fam = {f["query"]: (f["raw_occurrences"], f["length"]) for f in load(a.old, "all_native_families.json")}
    new_fam = {f["query"]: (f["raw_occurrences"], f["length"]) for f in load(a.new, "all_native_families.json")}
    if old_fam != new_fam:
        diffs.append(dict(kind="native_families", only_old=sorted(set(old_fam) - set(new_fam)),
                          only_new=sorted(set(new_fam) - set(old_fam)),
                          changed=sorted(k for k in set(old_fam) & set(new_fam) if old_fam[k] != new_fam[k])))
    old_m = {r["member"]: r for r in load(a.old, "member_boundaries.json")}
    new_m = {r["member"]: r for r in load(a.new, "member_boundaries.json")}
    for m in sorted(set(old_m) | set(new_m)):
        if m not in old_m or m not in new_m:
            diffs.append(dict(kind="member_presence", member=m, in_old=m in old_m, in_new=m in new_m))
            continue
        changed = {f: (old_m[m].get(f), new_m[m].get(f)) for f in FIELDS if old_m[m].get(f) != new_m[m].get(f)}
        if changed:
            diffs.append(dict(kind="member_fields", member=m, changed=changed))
    for name in FASTAS:
        o, n = fasta(a.old / "sequences" / name), fasta(a.new / "sequences" / name)
        bad = sorted(k for k in set(o) | set(n) if o.get(k) != n.get(k))
        if bad:
            diffs.append(dict(kind="sequence", file=name, records=bad[:50], count=len(bad)))
    report = dict(old=str(a.old), new=str(a.new), native_families=len(new_fam), members_old=len(old_m),
                  members_new=len(new_m), fields_compared=FIELDS, sequence_files_compared=FASTAS,
                  identical=not diffs, differences=diffs)
    a.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("native_families", "members_old", "members_new", "identical")}))
    sys.exit(0 if not diffs else 1)


if __name__ == "__main__":
    main()
