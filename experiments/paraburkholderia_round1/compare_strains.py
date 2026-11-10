"""Cross-strain X families and host positions for the Paraburkholderia round-1 panel.

Inputs per strain: a multisite-insertion run (member_boundaries.json and
sequences/candidate_insertion.fasta), the read_callout evidence table and the
assembly FASTA. Steps follow docs/PLAN_PARABURKHOLDERIA_MULTI_ISOLATE.md:

P2  strict X families: blastn all-vs-all of candidate X; a pair is strict when
    one HSP covers >=90% of both sequences at >=95% identity. A connected
    component is a strict group only if every pair is strict; otherwise it is
    clustering_ambiguous and its members stay in singleton groups for counting.
    Within a strain, members of one strict group whose candidate intervals
    overlap reciprocally >=90% are one locus (all pairs required).
P3  host positions: each locus's genomic 500 bp flanks are placed in the other
    strain through nucmer --maxmatch blocks. A flank is placed only by blocks
    covering >=80% of it and only if all such blocks agree on one location.
    Unplaced flanks make the locus unresolved; unresolved is never counted as
    a new position.

Read support uses the v1 read_callout classes (full-span molecules per member);
the audited P4 read rule is not implemented here and is reported as such.
"""

import argparse
import csv
import json
import subprocess
from collections import defaultdict
from pathlib import Path

FAMILY_COVER = 0.90
FAMILY_IDENTITY = 95.0
LOCUS_OVERLAP = 0.90
FLANK_BP = 500
MIN_FLANK_BP = 250
FLANK_BLOCK_COVER = 0.80
SAME_PLACE_BP = 1000  # placements of one flank within this distance are one location
EMPTY_GAP_BP = 50  # projected flank ends this close in the other strain: no insert there
SHARED_SLOP_BP = 100
MAX_SPAN_BP = 50000  # projected flanks farther apart than this: not one position
MIN_SUPPORT = 2  # v1 full-span molecules for a locus to count as read supported


def run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def load_strain(strain, run_dir, evidence_dir, fasta):
    rows = json.loads((run_dir / "results/member_boundaries.json").read_text())
    xseq = read_fasta(run_dir / "sequences/candidate_insertion.fasta")
    support = {}
    if evidence_dir and (evidence_dir / "member_read_support.tsv").exists():
        for r in csv.DictReader(open(evidence_dir / "member_read_support.tsv"), delimiter="\t"):
            support[r["member"]] = (int(r["full"]), r["class"])
    lengths = {k: len(v) for k, v in read_fasta(fasta).items()}
    members = []
    for r in rows:
        if r.get("candidate_start0") is None or r["member"] not in xseq:
            continue
        s, e = r["candidate_start0"], r["candidate_end0"]
        full, cls = support.get(r["member"], (0, "no_read_data"))
        members.append(dict(
            strain=strain, member=r["member"], family=r["family"], method=r["method"],
            contig=r["accession"], s=s, e=e, length=e - s, strand=r["strand"],
            boundary_source=r["boundary_source"],
            up_len=min(FLANK_BP, s), down_len=min(FLANK_BP, lengths[r["accession"]] - e),
            full=full, read_class=cls, uid=strain + "|" + r["member"],
        ))
    return members, xseq, lengths


def strict_family_groups(members, xseqs, work):
    fa = work / "all_x.fasta"
    with open(fa, "w") as out:
        for m in members:
            out.write(">%s\n%s\n" % (m["uid"], xseqs[m["strain"]][m["member"]]))
    hits = work / "x_vs_x.tsv"
    run(["blastn", "-task", "blastn", "-dust", "no", "-query", str(fa), "-subject", str(fa),
         "-evalue", "1e-10", "-max_target_seqs", "100000", "-max_hsps", "1000",
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send qlen slen", "-out", str(hits)])
    strict = defaultdict(set)
    pairs = []
    for line in open(hits):
        q, t, pid, alen, qs, qe, ss, se, ql, sl = line.split("\t")
        if q == t:
            continue
        qcov = (int(qe) - int(qs) + 1) / int(ql)
        scov = (abs(int(se) - int(ss)) + 1) / int(sl)
        if float(pid) >= FAMILY_IDENTITY and qcov >= FAMILY_COVER and scov >= FAMILY_COVER:
            strict[q].add(t)
            strict[t].add(q)
            pairs.append((q, t, float(pid), round(qcov, 3), round(scov, 3)))
    # Connected components; keep only cliques as strict groups.
    seen, groups, ambiguous = set(), [], []
    for m in members:
        u = m["uid"]
        if u in seen:
            continue
        comp, stack = set(), [u]
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(strict[x] - comp)
        seen |= comp
        clique = all(b in strict[a] for a in comp for b in comp if a != b)
        if clique:
            groups.append(sorted(comp))
        else:
            ambiguous.append(sorted(comp))
            groups.extend([x] for x in sorted(comp))
    return groups, ambiguous, pairs


def collapse_loci(members, groups):
    by_uid = {m["uid"]: m for m in members}
    loci, locus_ambiguous = [], 0
    for gi, grp in enumerate(groups):
        per_strain = defaultdict(list)
        for u in grp:
            per_strain[by_uid[u]["strain"]].append(by_uid[u])
        for strain, ms in per_strain.items():
            def same(a, b):
                ov = min(a["e"], b["e"]) - max(a["s"], b["s"])
                return a["contig"] == b["contig"] and ov >= LOCUS_OVERLAP * max(a["length"], b["length"])
            remaining = sorted(ms, key=lambda m: (m["contig"], m["s"]))
            while remaining:
                seed = remaining.pop(0)
                cluster = [seed] + [m for m in remaining if same(seed, m)]
                if not all(same(a, b) for a in cluster for b in cluster if a is not b):
                    locus_ambiguous += 1
                    cluster = [seed]
                remaining = [m for m in remaining if m not in cluster]
                rep = max(cluster, key=lambda m: (m["boundary_source"] != "native_interval_only_unresolved",
                                                  m["up_len"] + m["down_len"], -abs(m["length"] - sorted(x["length"] for x in cluster)[len(cluster) // 2]), m["member"]))
                loci.append(dict(
                    group=gi, strain=strain, contig=rep["contig"], s=rep["s"], e=rep["e"], length=rep["length"],
                    members=[m["member"] for m in cluster], rep=rep["member"],
                    flanks_complete=rep["up_len"] >= MIN_FLANK_BP and rep["down_len"] >= MIN_FLANK_BP,
                    full=max(m["full"] for m in cluster), read_classes=sorted({m["read_class"] for m in cluster}),
                ))
    for i, l in enumerate(loci):
        l["locus_id"] = "%s:L%d" % (l["strain"], i)
        l["supported"] = l["full"] >= MIN_SUPPORT
    return loci, locus_ambiguous


def nucmer_blocks(ref_fa, qry_fa, prefix):
    run(["nucmer", "--maxmatch", "-p", str(prefix), str(ref_fa), str(qry_fa)], stderr=subprocess.DEVNULL)
    out = subprocess.run(["show-coords", "-rclTH", str(prefix) + ".delta"], check=True, capture_output=True, text=True).stdout
    blocks = defaultdict(list)
    for line in out.splitlines():
        f = line.split("\t")
        s1, e1, s2, e2 = map(int, f[:4])
        blocks[f[11]].append((s1 - 1, e1, s2, e2, f[12]))  # ref 0-based half-open; query 1-based, may be reversed
    for v in blocks.values():
        v.sort()
    return blocks


def place(blocks, contig, start, end, inner):
    """Project `inner` (a ref coordinate at the flank's X-facing end) into the other strain."""
    need = FLANK_BLOCK_COVER * (end - start)
    hits = []
    for s1, e1, s2, e2, qname in blocks.get(contig, ()):
        if s1 >= end:
            break
        if min(e1, end) - max(s1, start) >= need:
            fwd = s2 <= e2
            pos = s2 - 1 + (inner - s1) if fwd else s2 - (inner - s1)
            hits.append((qname, pos, fwd))
    places = []
    for h in hits:
        if not any(h[0] == p[0] and abs(h[1] - p[1]) <= SAME_PLACE_BP and h[2] == p[2] for p in places):
            places.append(h)
    if not places:
        return None, "unplaced"
    if len(places) > 1:
        return None, "multiple_places"
    return places[0], None


def compare_pair(a, b, loci, blocks_ab):
    """Classify every locus of strain a by what strain b has at the same host position."""
    b_loci = defaultdict(list)
    for l in loci:
        if l["strain"] == b:
            b_loci[(l["group"], l["contig"])].append(l)
    out = []
    for l in (x for x in loci if x["strain"] == a):
        rec = dict(locus_id=l["locus_id"], other=b, status=None, reason="", other_contig="", other_start0="", other_end0="", matched_locus="")
        if not l["flanks_complete"]:
            rec.update(status="unresolved", reason="flank_short_or_contig_edge")
            out.append(rec)
            continue
        up, why_up = place(blocks_ab, l["contig"], l["s"] - FLANK_BP, l["s"], l["s"])
        down, why_down = place(blocks_ab, l["contig"], l["e"], l["e"] + FLANK_BP, l["e"])
        if up is None or down is None:
            rec.update(status="unresolved", reason="up_" + why_up if up is None else "down_" + why_down)
            out.append(rec)
            continue
        if up[0] != down[0] or up[2] != down[2]:
            rec.update(status="unresolved", reason="flanks_on_different_contigs_or_strands")
            out.append(rec)
            continue
        lo, hi = sorted((up[1], down[1]))
        if (up[1] > down[1]) == up[2] or hi - lo > MAX_SPAN_BP:
            rec.update(status="unresolved", reason="flank_order_or_distance_inconsistent")
            out.append(rec)
            continue
        rec.update(other_contig=up[0], other_start0=lo, other_end0=hi)
        match = [m for m in b_loci.get((l["group"], up[0]), ()) if m["s"] < hi + SHARED_SLOP_BP and m["e"] > lo - SHARED_SLOP_BP]
        if match:
            rec.update(status="shared_position", matched_locus=",".join(m["locus_id"] for m in match))
        elif hi - lo <= EMPTY_GAP_BP:
            rec.update(status="empty_in_other_assembly")
        else:
            rec.update(status="other_sequence_in_other_assembly", reason="gap_%dbp" % (hi - lo))
        out.append(rec)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strain", nargs=4, action="append", required=True, metavar=("NAME", "RUN_DIR", "EVIDENCE_DIR", "FASTA"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    work = args.out / "work"
    work.mkdir()

    members, xseqs, fastas = [], {}, {}
    for name, run_dir, ev_dir, fa in args.strain:
        ms, xs, _ = load_strain(name, Path(run_dir), Path(ev_dir) if ev_dir != "-" else None, Path(fa))
        members += ms
        xseqs[name] = xs
        fastas[name] = Path(fa)
    groups, ambiguous, pairs = strict_family_groups(members, xseqs, work)
    loci, locus_ambiguous = collapse_loci(members, groups)
    strains = [s[0] for s in args.strain]

    comparisons = []
    for a in strains:
        for b in strains:
            if a != b:
                blocks = nucmer_blocks(fastas[a], fastas[b], work / ("%s_vs_%s" % (a, b)))
                comparisons += [dict(strain=a, **r) for r in compare_pair(a, b, loci, blocks)]

    by_id = {l["locus_id"]: l for l in loci}
    with open(args.out / "loci.tsv", "w") as h:
        cols = ["locus_id", "strain", "group", "contig", "s", "e", "length", "rep", "members", "flanks_complete", "full", "supported", "read_classes"]
        h.write("\t".join(cols) + "\n")
        for l in loci:
            h.write("\t".join(",".join(map(str, l[c])) if isinstance(l[c], list) else str(l[c]) for c in cols) + "\n")
    with open(args.out / "position_comparisons.tsv", "w") as h:
        cols = ["strain", "locus_id", "other", "status", "reason", "other_contig", "other_start0", "other_end0", "matched_locus"]
        h.write("\t".join(cols) + "\n")
        for r in comparisons:
            h.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(args.out / "family_pair_matches.tsv", "w") as h:
        h.write("a\tb\tidentity\tcover_a\tcover_b\n")
        for p in pairs:
            h.write("\t".join(map(str, p)) + "\n")

    # Per strict group shared by >=2 strains: what does adding strain b give relative to strain a?
    summary = []
    for gi, grp in enumerate(groups):
        present = sorted({by_id[l["locus_id"]]["strain"] for l in loci if l["group"] == gi})
        if len(present) < 2:
            continue
        lengths = sorted(l["length"] for l in loci if l["group"] == gi)
        for b in present:
            for a in present:
                if a == b:
                    continue
                rows = [r for r in comparisons if r["strain"] == b and r["other"] == a and by_id[r["locus_id"]]["group"] == gi]
                st = lambda s: [r for r in rows if r["status"] == s]
                new = [r for r in rows if r["status"] in ("empty_in_other_assembly", "other_sequence_in_other_assembly")]
                summary.append(dict(
                    group=gi, x_length_min=lengths[0], x_length_max=lengths[-1], added_strain=b, baseline_strain=a,
                    added_loci=len(rows), added_loci_read_supported=sum(by_id[r["locus_id"]]["supported"] for r in rows),
                    shared_position=len(st("shared_position")),
                    new_position_empty_in_baseline=len(st("empty_in_other_assembly")),
                    new_position_other_sequence_in_baseline=len(st("other_sequence_in_other_assembly")),
                    new_position_read_supported=sum(by_id[r["locus_id"]]["supported"] for r in new),
                    unresolved=len(st("unresolved")),
                ))
    with open(args.out / "family_increment.tsv", "w") as h:
        cols = list(summary[0]) if summary else ["group"]
        h.write("\t".join(cols) + "\n")
        for r in summary:
            h.write("\t".join(str(r[c]) for c in cols) + "\n")
    info = dict(
        strains=strains, members=len(members), strict_groups=len(groups) - sum(len(c) for c in ambiguous),
        clustering_ambiguous_components=len(ambiguous), loci=len(loci), locus_overlap_ambiguous=locus_ambiguous,
        groups_in_two_or_more_strains=len({r["group"] for r in summary}),
        read_rule="v1 read_callout full-span molecules >= %d per locus (audited P4 rule not implemented)" % MIN_SUPPORT,
        thresholds=dict(FAMILY_COVER=FAMILY_COVER, FAMILY_IDENTITY=FAMILY_IDENTITY, LOCUS_OVERLAP=LOCUS_OVERLAP,
                        FLANK_BP=FLANK_BP, MIN_FLANK_BP=MIN_FLANK_BP, FLANK_BLOCK_COVER=FLANK_BLOCK_COVER,
                        SAME_PLACE_BP=SAME_PLACE_BP, EMPTY_GAP_BP=EMPTY_GAP_BP, SHARED_SLOP_BP=SHARED_SLOP_BP,
                        MAX_SPAN_BP=MAX_SPAN_BP, MIN_SUPPORT=MIN_SUPPORT),
    )
    (args.out / "summary.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
