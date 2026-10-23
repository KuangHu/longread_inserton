"""Audited filled-site read check for assembly insertion candidates (protocol section 4).

Reads of one sample are mapped to that sample's complete assembly with minimap2
(-c, secondaries kept). For every exported candidate X this script decides,
read by read, whether one continuous alignment spans X plus 100 bp anchors on
both sides with enough aligned reference bases, and whether that placement is
contested by another alignment of the same read. It never changes discovery
output and never upgrades boundary status.

Subcommands:
  cap-reads  list reads touching candidates whose secondary count hit the -N cap
  extract    write those reads from the FASTQ for a sensitivity re-mapping
  classify   per-read evidence, per-member audited class, v1 differences, summary

Output schema: audited_filled_site_v1. It is not the v1 read_callout format.
"""

import argparse
import bisect
import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCHEMA = "audited_filled_site_v1"
ANCHOR_BP = 100
MIN_ANCHOR_COVER = 0.90  # aligned (M/=/X) reference bases in each 100 bp anchor
MIN_X_COVER = 0.90  # aligned reference bases inside the candidate
MAX_X_INDEL = 0.10  # I+D bases inside the candidate, as a fraction of its length
MIN_UNITS = 2  # identity-distinguishable support units for audited full span
X_DELETED = 0.50  # one deletion of >= this fraction of X inside X: read lacks X
COMPETE_AS = 0.95  # another alignment of the same read region with AS >= this ratio competes
COMPETE_QUERY_OVERLAP = 0.50
SECONDARY_CAP = 50  # minimap2 -N used for the primary mapping
RECHECK_CAP = 200  # -N used for the sensitivity re-mapping
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

CLASSES = [
    "edge_limited",
    "audited_full_span",
    "audited_full_span_with_conflict",
    "single_unit_full_span",
    "ambiguous_competing",
    "full_span_cigar_fail",
    "conflict_x_absent",
    "one_side_only",
    "internal_only",
    "no_coverage",
]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_paf(line):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        raise ValueError("Malformed PAF line: " + line[:80])
    tags = {t[:2]: t[5:] for t in f[12:]}
    return dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], tlen=int(f[6]),
                ts=int(f[7]), te=int(f[8]), mapq=int(f[11]), tp=tags.get("tp", "P"),
                AS=int(tags.get("AS", 0)), cg=tags.get("cg", ""))


def fasta_lengths(path):
    lengths, name = {}, None
    with open_text(path) as h:
        for line in h:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def load_members(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as h:
        rows = json.load(h)
    out = []
    for r in rows:
        if r.get("candidate_start0") is None:
            continue
        out.append(dict(member=r["member"], family=r["family"], method=r["method"], contig=r["accession"],
                        s=r["candidate_start0"], e=r["candidate_end0"], boundary_source=r["boundary_source"]))
    return out


class Windows:
    """Candidate windows (X plus anchors) indexed by contig for overlap queries."""

    def __init__(self, members):
        self.by = defaultdict(list)
        for m in members:
            self.by[m["contig"]].append((m["s"] - ANCHOR_BP, m["e"] + ANCHOR_BP, m))
        for v in self.by.values():
            v.sort(key=lambda x: x[0])
        self.starts = {c: [w[0] for w in v] for c, v in self.by.items()}
        self.span = max((w[1] - w[0] for v in self.by.values() for w in v), default=0)

    def hits(self, contig, ts, te):
        v = self.by.get(contig, [])
        i = bisect.bisect_left(self.starts.get(contig, []), ts - self.span)
        out = []
        for ws, we, m in v[i:]:
            if ws >= te:
                break
            if we > ts:
                out.append(m)
        return out


def cigar_metrics(a, s, e):
    """Aligned/indel bases of one alignment over [s-100,s), [s,e), [e,e+100)."""
    regions = dict(left=(s - ANCHOR_BP, s), x=(s, e), right=(e, e + ANCHOR_BP))
    cov = dict.fromkeys(regions, 0)
    ins_x = del_x = boundary_ins = n_window = largest_del_x = 0
    pos = a["ts"]
    for length, op in CIGAR.findall(a["cg"]):
        n = int(length)
        if op in "M=X":
            for k, (rs, re_) in regions.items():
                cov[k] += max(0, min(pos + n, re_) - max(pos, rs))
            pos += n
        elif op in "DN":
            ov = max(0, min(pos + n, e) - max(pos, s))
            if op == "N":
                n_window += max(0, min(pos + n, e + ANCHOR_BP) - max(pos, s - ANCHOR_BP))
            else:
                del_x += ov
                largest_del_x = max(largest_del_x, ov)
            pos += n
        elif op == "I":
            # An insertion sits between reference pos-1 and pos: inside X only if s < pos < e.
            if s < pos < e:
                ins_x += n
            elif pos in (s, e):
                boundary_ins = max(boundary_ins, n)
    return dict(cov_left=cov["left"], cov_x=cov["x"], cov_right=cov["right"], ins_x=ins_x, del_x=del_x,
                largest_del_x=largest_del_x, boundary_ins=boundary_ins, cigar_N_in_window=n_window)


def read_alignments(paf, windows):
    """Two passes: reads touching a candidate window, then every alignment of those reads."""
    touching = set()
    with open_text(paf) as h:
        for line in h:
            a = parse_paf(line)
            if windows.hits(a["t"], a["ts"], a["te"]):
                touching.add(a["q"])
    by_read = defaultdict(list)
    with open_text(paf) as h:
        for line in h:
            q = line.split("\t", 1)[0]
            if q in touching:
                by_read[q].append(parse_paf(line))
    return by_read


def competing(a, alns, L):
    """Another alignment of the same read region elsewhere with a near-equal score."""
    for b in alns:
        if b is a:
            continue
        same_place = b["t"] == a["t"] and abs(b["ts"] - a["ts"]) <= L + 200 and b["strand"] == a["strand"]
        if same_place:
            continue
        ov = min(a["qe"], b["qe"]) - max(a["qs"], b["qs"])
        if ov >= COMPETE_QUERY_OVERLAP * (a["qe"] - a["qs"]) and b["AS"] >= COMPETE_AS * a["AS"]:
            return "%s:%d-%d(%s,AS=%d)" % (b["t"], b["ts"], b["te"], b["tp"], b["AS"])
    return ""


def classify_read(m, alns, contig_len, cap_reads, recheck):
    """Return per-alignment rows and the read's verdict for one candidate."""
    s, e = m["s"], m["e"]
    L = e - s
    pool = recheck.get(alns[0]["q"], alns) if alns and alns[0]["q"] in cap_reads else alns
    n_secondary = sum(1 for x in pool if x["tp"] != "P")
    capped = alns[0]["q"] in cap_reads and alns[0]["q"] not in recheck
    recheck_capped = alns[0]["q"] in recheck and n_secondary >= RECHECK_CAP
    rows, verdicts = [], []
    for a in pool:
        if a["t"] != m["contig"] or a["te"] <= s - ANCHOR_BP or a["ts"] >= e + ANCHOR_BP:
            continue
        met = cigar_metrics(a, s, e)
        spans = a["ts"] <= s - ANCHOR_BP and a["te"] >= e + ANCHOR_BP
        reasons = []
        if met["largest_del_x"] >= X_DELETED * L:
            # A read lacking X is conflict evidence only from its best placement; a secondary
            # alignment over a neighbouring copy is a competing placement, not an absent X.
            if a["tp"] == "P":
                status = "x_deleted"
            else:
                status, reasons = "competing", ["secondary_alignment_lacks_x"]
        elif spans:
            if met["cov_left"] < MIN_ANCHOR_COVER * ANCHOR_BP:
                reasons.append("left_anchor_cover")
            if met["cov_right"] < MIN_ANCHOR_COVER * ANCHOR_BP:
                reasons.append("right_anchor_cover")
            if met["cov_x"] < MIN_X_COVER * L:
                reasons.append("x_cover")
            if met["ins_x"] + met["del_x"] > MAX_X_INDEL * L:
                reasons.append("x_indel")
            if met["cigar_N_in_window"]:
                reasons.append("cigar_N")
            if reasons:
                status = "full_span_cigar_fail"
            elif a["tp"] != "P":
                status, reasons = "competing", ["secondary_at_candidate"]
            else:
                other = competing(a, pool, L)
                if other:
                    status, reasons = "competing", ["near_equal_alignment " + other]
                elif capped or recheck_capped:
                    status, reasons = "competing", ["secondary_cap_unresolved"]
                else:
                    status = "full_support"
        elif a["ts"] <= s - ANCHOR_BP and met["cov_left"] >= MIN_ANCHOR_COVER * ANCHOR_BP:
            status = "left_only"
        elif a["te"] >= e + ANCHOR_BP and met["cov_right"] >= MIN_ANCHOR_COVER * ANCHOR_BP:
            status = "right_only"
        else:
            status = "internal"
        verdicts.append(status)
        rows.append(dict(read=a["q"], tp=a["tp"], mapq=a["mapq"], AS=a["AS"], strand=a["strand"], qs=a["qs"], qe=a["qe"],
                         ts=a["ts"], te=a["te"], n_secondary=n_secondary, cap_reached=capped or recheck_capped,
                         rechecked=a["q"] in recheck, status=status, reason=";".join(reasons), **met))
    order = ["full_support", "x_deleted", "competing", "full_span_cigar_fail", "left_only", "right_only", "internal"]
    verdict = min(verdicts, key=order.index) if verdicts else None
    return rows, verdict


def unit_key(rows):
    """Reads with an identical alignment signature at this locus are one support unit."""
    best = min(rows, key=lambda r: (r["status"] != "full_support", r["ts"]))
    return (best["ts"], best["te"], best["strand"], best["qe"] - best["qs"], best["cov_x"], best["ins_x"], best["del_x"])


def classify_member(m, unit_verdicts, contig_len):
    if m["s"] - ANCHOR_BP < 0 or m["e"] + ANCHOR_BP > contig_len.get(m["contig"], 0):
        return "edge_limited"
    n = Counter(unit_verdicts.values())
    if n["full_support"] >= MIN_UNITS:
        return "audited_full_span_with_conflict" if n["x_deleted"] >= MIN_UNITS else "audited_full_span"
    if n["full_support"] == 1:
        return "single_unit_full_span"
    if n["competing"]:
        return "ambiguous_competing"
    if n["full_span_cigar_fail"]:
        return "full_span_cigar_fail"
    if n["x_deleted"]:
        return "conflict_x_absent"
    if n["left_only"] or n["right_only"]:
        return "one_side_only"
    if n["internal"]:
        return "internal_only"
    return "no_coverage"


def cap_reads(by_read):
    return {q for q, alns in by_read.items() if sum(1 for a in alns if a["tp"] != "P") >= SECONDARY_CAP}


def classify(args):
    members = load_members(args.members)
    contig_len = fasta_lengths(args.reference)
    windows = Windows(members)
    by_read = read_alignments(args.paf, windows)
    capped = cap_reads(by_read)
    recheck = defaultdict(list)
    if args.recheck_paf:
        for line in open_text(args.recheck_paf):
            a = parse_paf(line)
            recheck[a["q"]].append(a)
    per_member = defaultdict(dict)  # member -> read -> rows
    verdicts = defaultdict(dict)
    by_id = {m["member"]: m for m in members}
    for q, alns in by_read.items():
        targets = {m["member"]: m for a in alns for m in windows.hits(a["t"], a["ts"], a["te"])}
        for mid, m in targets.items():
            rows, verdict = classify_read(m, alns, contig_len, capped, recheck)
            if verdict:
                per_member[mid][q] = rows
                verdicts[mid][q] = verdict

    args.out.mkdir(parents=True, exist_ok=False)
    v1 = {}
    if args.v1 and Path(args.v1).exists():
        v1 = {r["member"]: r for r in csv.DictReader(open(args.v1), delimiter="\t")}
    cols = ["member", "read", "unit", "duplicate_of_unit", "tp", "mapq", "AS", "strand", "qs", "qe", "ts", "te",
            "cov_left", "cov_x", "cov_right", "ins_x", "del_x", "largest_del_x", "boundary_ins", "cigar_N_in_window",
            "n_secondary", "cap_reached", "rechecked", "status", "reason"]
    summary_rows = []
    with open(args.out / "read_evidence.tsv", "w") as ev, open(args.out / "member_audit.tsv", "w") as mo:
        ev.write("\t".join(cols) + "\n")
        mcols = ["member", "family", "method", "contig", "start0", "end0", "length", "boundary_source",
                 "units_full_support", "units_x_deleted", "units_competing", "units_cigar_fail", "units_left_only",
                 "units_right_only", "units_internal", "reads_total", "duplicate_reads", "cap_reads",
                 "audited_class", "independence", "assembly_consistency"]
        mo.write("\t".join(mcols) + "\n")
        for m in members:
            reads = per_member.get(m["member"], {})
            units, unit_verdict, dups = {}, {}, 0
            for q in sorted(reads):
                key = unit_key(reads[q])
                if key in units:
                    dups += 1
                    dup_of = units[key]
                else:
                    units[key] = q
                    dup_of = ""
                    unit_verdict[q] = verdicts[m["member"]][q]
                for r in reads[q]:
                    ev.write("\t".join(str(x) for x in [m["member"], q, units[key], dup_of] + [r[c] for c in cols[4:]]) + "\n")
            cls = classify_member(m, unit_verdict, contig_len)
            n = Counter(unit_verdict.values())
            row = dict(member=m["member"], family=m["family"], method=m["method"], contig=m["contig"], start0=m["s"],
                       end0=m["e"], length=m["e"] - m["s"], boundary_source=m["boundary_source"],
                       units_full_support=n["full_support"], units_x_deleted=n["x_deleted"], units_competing=n["competing"],
                       units_cigar_fail=n["full_span_cigar_fail"], units_left_only=n["left_only"],
                       units_right_only=n["right_only"], units_internal=n["internal"], reads_total=len(reads),
                       duplicate_reads=dups, cap_reads=sum(1 for q in reads if q in capped),
                       audited_class=cls, independence=args.independence, assembly_consistency=args.assembly_consistency)
            summary_rows.append(row)
            mo.write("\t".join(str(row[c]) for c in mcols) + "\n")
    with open(args.out / "v1_vs_audited.tsv", "w") as h:
        h.write("member\tv1_class\tv1_full\taudited_class\taudited_full_units\tchanged\n")
        for r in summary_rows:
            old = v1.get(r["member"], {})
            v1_ok = old.get("class") in ("full_span", "full_span_and_empty_allele")
            new_ok = r["audited_class"].startswith("audited_full_span")
            h.write("\t".join(map(str, [r["member"], old.get("class", "not_in_v1"), old.get("full", ""), r["audited_class"],
                                         r["units_full_support"], v1_ok != new_ok])) + "\n")
    counts = Counter(r["audited_class"] for r in summary_rows)
    long_rows = [r for r in summary_rows if r["length"] >= 400]
    summary = dict(
        schema=SCHEMA, sample=args.sample, members=len(summary_rows), reads_touching_candidates=len(by_read),
        cap_reads=len(capped), cap_reads_rechecked=len(set(recheck) & capped),
        audited_class_counts={c: counts.get(c, 0) for c in CLASSES},
        long_member_class_counts={c: sum(1 for r in long_rows if r["audited_class"] == c) for c in CLASSES},
        v1_comparison=dict(v1_members=len(v1), changed=sum(
            (v1.get(r["member"], {}).get("class") in ("full_span", "full_span_and_empty_allele")) != r["audited_class"].startswith("audited_full_span")
            for r in summary_rows)),
        independence=args.independence, assembly_consistency=args.assembly_consistency,
        thresholds=dict(ANCHOR_BP=ANCHOR_BP, MIN_ANCHOR_COVER=MIN_ANCHOR_COVER, MIN_X_COVER=MIN_X_COVER,
                        MAX_X_INDEL=MAX_X_INDEL, MIN_UNITS=MIN_UNITS, X_DELETED=X_DELETED, COMPETE_AS=COMPETE_AS,
                        COMPETE_QUERY_OVERLAP=COMPETE_QUERY_OVERLAP, SECONDARY_CAP=SECONDARY_CAP, RECHECK_CAP=RECHECK_CAP),
        notes=["Boundary status is never upgraded by read support.",
               "Support units are distinct reads after collapsing identical alignment signatures; "
               "they are not certified independent molecules unless independence says so."],
    )
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("sample", "members", "cap_reads", "audited_class_counts")}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cap-reads")
    c.add_argument("--members", type=Path, required=True)
    c.add_argument("--paf", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    x = sub.add_parser("extract")
    x.add_argument("--reads", type=Path, required=True)
    x.add_argument("--fastq", type=Path, required=True)
    x.add_argument("--out", type=Path, required=True)
    k = sub.add_parser("classify")
    k.add_argument("--members", type=Path, required=True, help="member_boundaries.json(.gz)")
    k.add_argument("--reference", type=Path, required=True)
    k.add_argument("--paf", type=Path, required=True)
    k.add_argument("--recheck-paf", type=Path)
    k.add_argument("--v1", type=Path, help="v1 member_read_support.tsv for the difference table")
    k.add_argument("--sample", required=True)
    k.add_argument("--independence", default="independence_unverified")
    k.add_argument("--assembly-consistency", default="unknown")
    k.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    try:
        if a.cmd == "cap-reads":
            capped = cap_reads(read_alignments(a.paf, Windows(load_members(a.members))))
            a.out.write_text("".join(q + "\n" for q in sorted(capped)))
            print(len(capped))
        elif a.cmd == "extract":
            wanted = set(a.reads.read_text().split())
            with open_text(a.fastq) as h, open(a.out, "w") as out:
                while True:
                    rec = [h.readline() for _ in range(4)]
                    if not rec[0]:
                        break
                    if rec[0][1:].split()[0] in wanted:
                        out.writelines(rec)
        else:
            classify(a)
    except (ValueError, OSError, KeyError) as error:
        p.exit(1, "Error: %s\n" % error)


if __name__ == "__main__":
    main()
