"""Fixed-rule spot-check table for the five-sample acceptance run (protocol T1).

Selection (fixed, by member ID order): all AR20-38 PaAR2038_RS_R1 members; per
sample the first 5 priority candidates, then the first short (<400 bp), first
edge/censored, first tnp-unsupported and first conflict/ambiguous member if
present. For each: re-slice X and both flanks from the assembly and compare with
the exports, boundary evidence, and audited read evidence. The verdict column
is filled by the reviewer after reading the row.
"""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CONFLICT = ("conflict_x_absent", "audited_full_span_with_conflict", "ambiguous_competing")


def fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def rc(s):
    return s.translate(COMP)[::-1]


def select(rows):
    by_id = lambda r: (r["family"], int(r["member"].rsplit("_m", 1)[1]))
    rows = sorted(rows, key=by_id)
    picks = []
    if rows and rows[0]["sample"] == "PaAR2038":
        picks += [(r, "AR20-38 IS5 regression") for r in rows if r["family"] == "PaAR2038_RS_R1"]
    picks += [(r, "priority %d" % (i + 1)) for i, r in enumerate([r for r in rows if r["priority"] == "True"][:5])]
    rules = [("first short", lambda r: r["length_layer"] == "<400"),
             ("first edge/censored", lambda r: r["boundary_source"] == "homology_extent_censored_at_window"
              or r["audited_class"] == "edge_limited"),
             ("first tnp-unsupported", lambda r: r["tnp"] == "unsupported"),
             ("first conflict/ambiguous", lambda r: r["audited_class"] in CONFLICT)]
    for label, rule in rules:
        hit = next((r for r in rows if rule(r)), None)
        picks.append((hit, label) if hit else (None, label + ": none observed"))
    return picks


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True, help="sample_report.py output directory")
    p.add_argument("--proj", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout"))
    p.add_argument("--audit-root", default="callout_scale_v1")
    p.add_argument("--include-members", type=Path, help="file of extra member IDs to check (e.g. all class changes)")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    extra = set(a.include_members.read_text().split()) if a.include_members else set()
    rows = defaultdict(list)
    for r in csv.DictReader(open(a.report / "members.tsv"), delimiter="\t"):
        rows[r["sample"]].append(r)
    out = []
    for e in json.load(open(a.manifest))["samples"]:
        sample = e["sample_id"]
        genome = fasta(e["fasta"])
        seq_dir = Path(e["v1_callout"]) / "sequences"
        x, up, down = (fasta(seq_dir / f) for f in ("candidate_insertion.fasta", "upstream_500.fasta", "downstream_500.fasta"))
        bounds = {m["member"]: m for m in json.loads((Path(e["v1_callout"]) / "results/member_boundaries.json").read_text())}
        p_dir = a.proj / a.audit_root / sample / "evidence_audited"
        audit = {r["member"]: r for r in csv.DictReader(open(p_dir / "member_audit.tsv"), delimiter="\t")}
        reads = defaultdict(list)
        for r in csv.DictReader(open(p_dir / "read_evidence.tsv"), delimiter="\t"):
            reads[r["member"]].append(r)
        picks = select(rows[sample])
        chosen = {r["member"] for r, _ in picks if r}
        picks += [(r, "class changed in reclassification") for r in sorted(rows[sample], key=lambda r: r["member"])
                  if r["member"] in extra and r["member"] not in chosen]
        for r, why in picks:
            if r is None:
                out.append(dict(sample=sample, selection=why))
                continue
            m = bounds[r["member"]]
            g = genome[m["accession"]]
            s, t = m["candidate_start0"], m["candidate_end0"]
            core, left, right = g[s:t], g[max(0, s - 500):s], g[t:t + 500]
            if m["strand"] == "-":
                core, left, right = rc(core), rc(right), rc(left)
            slice_ok = (x.get(r["member"]) == core and up.get(r["member"]) == left and down.get(r["member"]) == right)
            ev = reads.get(r["member"], [])
            full = [q for q in ev if q["status"] == "full_support" and not q["duplicate_of_unit"]]
            best = max(full, key=lambda q: int(q["AS"]), default=None)
            competing = next((q["reason"] for q in ev if q["status"] == "competing"), "")
            a_row = audit[r["member"]]
            out.append(dict(
                sample=sample, selection=why, member=r["member"], family=r["family"], tnp=r["tnp"],
                slice="%s:%d-%d(%s) len %d" % (m["accession"], s, t, m["strand"], t - s),
                slice_matches_exports=slice_ok, flanks=f"{len(left)}/{len(right)}",
                boundary="%s; peers %s/%s; start spread %s; end spread %s; window touched %s" % (
                    m["boundary_source"], m.get("peer_count"), m.get("available_peers"),
                    m.get("oriented_start_spread_bp"), m.get("oriented_end_spread_bp"), m.get("context_window_touched")),
                v1="%s (full %s)" % (r["v1_class"], r["v1_full"]),
                audited="%s; units full %s, x_deleted %s, competing %s, cigar_fail %s; reads %s, duplicates %s, cap %s" % (
                    a_row["audited_class"], a_row["units_full_support"], a_row["units_x_deleted"], a_row["units_competing"],
                    a_row["units_cigar_fail"], a_row["reads_total"], a_row["duplicate_reads"], a_row["cap_reads"]),
                best_read="" if best is None else "cov_x %s/%d, ins %s del %s, anchors %s/%s, AS %s, MAPQ %s, %s" % (
                    best["cov_x"], t - s, best["ins_x"], best["del_x"], best["cov_left"], best["cov_right"],
                    best["AS"], best["mapq"], best["tp"]),
                competing_example=competing, verdict="",
            ))
    cols = ["sample", "selection", "member", "family", "tnp", "slice", "slice_matches_exports", "flanks", "boundary",
            "v1", "audited", "best_read", "competing_example", "verdict"]
    with open(a.out, "w") as h:
        h.write("\t".join(cols) + "\n")
        for r in out:
            h.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(len(out), "rows;", sum(1 for r in out if r.get("slice_matches_exports") is False), "slice mismatches")


if __name__ == "__main__":
    main()
