"""Classify long-read evidence for assembly-based insertion candidates.

Input: minimap2 PAF (with cg tags, secondaries kept) of molecule-named reads
against the complete matched assembly, plus a multisite-insertion run.

Evidence layers, all counted in distinct molecules (never subreads):

1. Per candidate member (reference loci from the assembly callout):
   full span, left/right junction, empty-allele reads (X absent), internal-only,
   ambiguous multimapping, no coverage, contig-edge limited.
2. Non-reference insertions, from two read patterns: a unique genomic anchor
   split-aligned to a copy of a candidate family away from every reference copy,
   or a large CIGAR insertion whose sequence (`insertions` + `extract`, then
   minimap2 back to the assembly) is a candidate-family copy. Inserted sequence
   matching no candidate is reported separately as unassigned. These are
   read-level observations, not validated transposition events.
3. Depth curves: layers 1-2 recomputed on deterministic molecule subsamples.

Subcommands: `insertions` lists large CIGAR insertions in unique alignments,
`extract` writes their sequences from the read FASTQ, `classify` combines all.

All thresholds are fixed constants below, chosen before looking at results.
"""

import argparse
import bisect
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

ANCHOR_BP = 100  # reference bases required on each side of a junction
MIN_MAPQ = 20  # unique placement for junction/anchor evidence
BIG_INDEL_FRACTION = 0.5  # of candidate length: deletion that removes X
EMPTY_GAP_BP = 200  # max read gap between flanks for split empty-allele reads
JUNCTION_SLOP_BP = 150  # breakpoint within this of a reference boundary = known site
CLUSTER_BP = 100  # merge non-reference breakpoints within this distance
ADJACENT_READ_GAP = 200  # anchor and X segments adjacent on the read
MIN_ANCHOR_LEN = 500
MIN_SITE_MOLECULES = 2
MIN_INSERTION_BP = 300  # CIGAR insertion length considered a candidate insert
MIN_INSERT_COVERED = 0.5  # fraction of inserted sequence aligned to a family copy
MIN_MEMBER_CONTAINED = 0.8  # or: fraction of a family copy contained in a larger insert
DEPTHS = (5, 10, 20, 40, 80)
SEEDS = (11, 29, 47)
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

CLASS_ORDER = [
    "edge_limited",
    "full_span_and_empty_allele",
    "empty_allele_only",
    "full_span",
    "both_junctions_separate",
    "one_junction",
    "ambiguous_multimapping",
    "internal_only",
    "no_coverage",
]


def molecule(qname):
    return qname.rsplit("/", 1)[0]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_paf(line):
    f = line.rstrip("\n").split("\t")
    tags = dict((t[:2], t[5:]) for t in f[12:])
    return dict(
        q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
        t=f[5], tlen=int(f[6]), ts=int(f[7]), te=int(f[8]),
        mapq=int(f[11]), tp=tags.get("tp", "P"), cg=tags.get("cg", ""),
    )


def indels(aln):
    """(ref_start, ref_end, op, length) for every I/D in the alignment."""
    pos = aln["ts"]
    out = []
    for length, op in CIGAR.findall(aln["cg"]):
        n = int(length)
        if op in "M=X":
            pos += n
        elif op == "D":
            out.append((pos, pos + n, "D", n))
            pos += n
        elif op == "I":
            out.append((pos, pos, "I", n))
    return out


def clean_span(aln, start, end, max_indel):
    """Alignment covers [start, end) with no indel >= max_indel inside it."""
    if aln["ts"] > start or aln["te"] < end:
        return False
    return not any(n >= max_indel and s < end and e >= start for s, e, _, n in indels(aln))


def unique(aln):
    return aln["tp"] == "P" and aln["mapq"] >= MIN_MAPQ


class MemberIndex:
    def __init__(self, members):
        self.by_contig = defaultdict(list)
        for m in members:
            self.by_contig[m["contig"]].append(m)
        for rows in self.by_contig.values():
            rows.sort(key=lambda m: m["s"])
        self.starts = {c: [m["s"] for m in rows] for c, rows in self.by_contig.items()}
        self.max_len = max((m["e"] - m["s"] for m in members), default=0)

    def overlapping(self, contig, start, end, pad=0):
        rows = self.by_contig.get(contig, [])
        lo = bisect.bisect_left(self.starts.get(contig, []), start - pad - self.max_len)
        out = []
        for m in rows[lo:]:
            if m["s"] - pad >= end:
                break
            if m["e"] + pad > start:
                out.append(m)
        return out


def load_members(run):
    rows = json.loads((run / "results/member_boundaries.json").read_text())
    families = json.loads((run / "results/all_families_annotated.json").read_text())
    tnp = {f["query"]: bool(f.get("tnp_supported")) for f in families}
    members = []
    for r in rows:
        if r.get("candidate_start0") is None:
            continue
        members.append(dict(
            member=r["member"], family=r["family"], method=r["method"],
            contig=r["accession"], s=r["candidate_start0"], e=r["candidate_end0"],
            boundary_source=r["boundary_source"], tnp_supported=tnp.get(r["family"], False),
        ))
    return members


def fasta_lengths(path):
    lengths, name = {}, None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            else:
                lengths[name] += len(line.strip())
    return lengths


def collect(paf, index):
    """Two passes: molecule lengths and reads touching candidates, then their alignments."""
    mol_len = {}
    touching = set()
    with open_text(paf) as handle:
        for line in handle:
            a = parse_paf(line)
            mol = molecule(a["q"])
            mol_len[mol] = max(mol_len.get(mol, 0), a["qlen"])
            if index.overlapping(a["t"], a["ts"], a["te"], pad=ANCHOR_BP + 2000):
                touching.add(a["q"])
    by_read = defaultdict(list)
    with open_text(paf) as handle:
        for line in handle:
            q = line.split("\t", 1)[0]
            if q in touching:
                by_read[q].append(parse_paf(line))
    return mol_len, by_read


def member_evidence(by_read, index, contig_len):
    """member -> molecule -> set of evidence flags."""
    ev = defaultdict(lambda: defaultdict(set))
    for q, alns in by_read.items():
        mol = molecule(q)
        for a in alns:
            for m in index.overlapping(a["t"], a["ts"], a["te"], pad=ANCHOR_BP):
                s, e = m["s"], m["e"]
                L = e - s
                big = max(50, int(BIG_INDEL_FRACTION * L))
                flags = ev[m["member"]][mol]
                u = unique(a)
                if clean_span(a, s - ANCHOR_BP, e + ANCHOR_BP, big):
                    flags.add("full" if u else "ambiguous")
                if clean_span(a, s - ANCHOR_BP, s + ANCHOR_BP, big):
                    flags.add("left" if u else "ambiguous")
                if clean_span(a, e - ANCHOR_BP, e + ANCHOR_BP, big):
                    flags.add("right" if u else "ambiguous")
                if u and a["ts"] <= s - ANCHOR_BP and a["te"] >= e + ANCHOR_BP:
                    for ds, de, op, n in indels(a):
                        ov = min(de, e) - max(ds, s)
                        if op == "D" and n >= BIG_INDEL_FRACTION * L and ov >= 0.5 * max(n, L):
                            flags.add("empty")
                if a["ts"] < e and a["te"] > s:
                    flags.add("internal")
        # Empty allele split across two unique alignments of one read piece.
        prim = [a for a in alns if unique(a)]
        for a in prim:
            for b in prim:
                if a is b or a["t"] != b["t"] or a["strand"] != b["strand"] or a["te"] > b["ts"] + EMPTY_GAP_BP:
                    continue
                gap = b["qs"] - a["qe"] if a["strand"] == "+" else a["qs"] - b["qe"]
                for m in index.overlapping(a["t"], a["te"], b["ts"]):
                    L = m["e"] - m["s"]
                    if (abs(a["te"] - m["s"]) <= JUNCTION_SLOP_BP and abs(b["ts"] - m["e"]) <= JUNCTION_SLOP_BP
                            and abs(gap) <= min(EMPTY_GAP_BP, 0.5 * L)):
                        ev[m["member"]][mol].add("empty")
    return ev


def classify(member, mols, contig_len):
    if member["s"] - ANCHOR_BP < 0 or member["e"] + ANCHOR_BP > contig_len.get(member["contig"], 0):
        return "edge_limited"
    n = lambda flag: sum(1 for f in mols.values() if flag in f)
    full, empty = n("full"), n("empty")
    if full and empty >= MIN_SITE_MOLECULES:
        return "full_span_and_empty_allele"
    if empty >= MIN_SITE_MOLECULES and not (full or n("left") or n("right")):
        return "empty_allele_only"
    if full:
        return "full_span"
    if n("left") and n("right"):
        return "both_junctions_separate"
    if n("left") or n("right"):
        return "one_junction"
    if n("ambiguous"):
        return "ambiguous_multimapping"
    if n("internal"):
        return "internal_only"
    return "no_coverage"


def nonreference_junctions(by_read, index):
    """Reads joining a unique anchor to a candidate-family copy away from known boundaries."""
    junctions = []
    known_splits = 0
    palindromic = 0
    for q, alns in by_read.items():
        prim = sorted((a for a in alns if a["tp"] == "P"), key=lambda a: a["qs"])
        if any(a["t"] == b["t"] and a["strand"] != b["strand"] and a["ts"] < b["te"] and b["ts"] < a["te"]
               for i, a in enumerate(prim) for b in prim[i + 1:]):
            palindromic += 1
        seg = []
        for a in prim:
            length = a["te"] - a["ts"]
            hits = index.overlapping(a["t"], a["ts"], a["te"])
            covered = max((min(a["te"], m["e"]) - max(a["ts"], m["s"]) for m in hits), default=0)
            if hits and covered >= 0.5 * length:
                fams = {m["family"] for m in hits if min(a["te"], m["e"]) - max(a["ts"], m["s"]) >= 0.5 * length}
                seg.append(("X", a, fams))
            elif unique(a) and length >= MIN_ANCHOR_LEN and covered < 0.2 * length:
                seg.append(("A", a, None))
            else:
                seg.append(("other", a, None))
        for (k1, u, f1), (k2, v, f2) in zip(seg, seg[1:]):
            if {k1, k2} != {"A", "X"} or abs(v["qs"] - u["qe"]) > ADJACENT_READ_GAP:
                continue
            anchor, fams = (u, f2) if k1 == "A" else (v, f1)
            anchor_first = k1 == "A"
            # Reference coordinate of the anchor end that touches X.
            at_end = anchor_first == (anchor["strand"] == "+")
            bp = anchor["te"] if at_end else anchor["ts"]
            near = [m for m in index.overlapping(anchor["t"], bp, bp + 1, pad=JUNCTION_SLOP_BP)
                    if m["family"] in fams and min(abs(bp - m["s"]), abs(bp - m["e"])) <= JUNCTION_SLOP_BP]
            if near:
                known_splits += 1
                continue
            junctions.append(dict(contig=anchor["t"], bp=bp, side="anchor_left" if at_end else "anchor_right",
                                  families=sorted(fams), molecule=molecule(q), read=q))
    return junctions, known_splits, palindromic


def scan_insertions(paf, out_path):
    """Large insertions in unique alignments, with read coordinates for extraction."""
    with open_text(paf) as handle, open(out_path, "w") as out:
        out.write("read\tqstart\tqend\tcontig\tbp\tlength\n")
        for line in handle:
            if "\tcg:Z:" not in line:
                continue
            a = parse_paf(line)
            if not unique(a):
                continue
            ref, rd = a["ts"], 0
            for length, op in CIGAR.findall(a["cg"]):
                n = int(length)
                if op == "I" and n >= MIN_INSERTION_BP:
                    q0 = a["qs"] + rd if a["strand"] == "+" else a["qe"] - rd - n
                    out.write("%s\t%d\t%d\t%s\t%d\t%d\n" % (a["q"], q0, q0 + n, a["t"], ref, n))
                if op in "M=XD":
                    ref += n
                if op in "M=XI":
                    rd += n


def extract_insertions(table, fastq, out_fa):
    wanted = defaultdict(list)
    with open(table) as handle:
        next(handle)
        for i, line in enumerate(handle):
            read, q0, q1 = line.split("\t")[:3]
            wanted[read].append((i, int(q0), int(q1)))
    with open_text(fastq) as handle, open(out_fa, "w") as out:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()
            handle.readline()
            name = header[1:].split()[0]
            for i, q0, q1 in wanted.get(name, ()):
                out.write(">ins%d %s:%d-%d\n%s\n" % (i, name, q0, q1, seq[q0:q1]))


def insertion_junctions(table, ins_paf, index):
    """Assign extracted insertions to candidate families via their assembly alignments."""
    rows = []
    with open(table) as handle:
        next(handle)
        for line in handle:
            read, q0, q1, contig, bp, length = line.rstrip("\n").split("\t")
            rows.append(dict(read=read, contig=contig, bp=int(bp), length=int(length)))
    fams = defaultdict(set)
    for line in open_text(ins_paf):
        a = parse_paf(line)
        if a["tp"] != "P":
            continue
        span = a["te"] - a["ts"]
        mostly_insert = (a["qe"] - a["qs"]) >= MIN_INSERT_COVERED * a["qlen"]
        for m in index.overlapping(a["t"], a["ts"], a["te"]):
            ov = min(a["te"], m["e"]) - max(a["ts"], m["s"])
            # Either the insert is mostly this copy, or it carries a near-full copy plus other sequence.
            if ov >= 0.5 * span and (mostly_insert or ov >= MIN_MEMBER_CONTAINED * (m["e"] - m["s"])):
                fams[int(a["q"][3:])].add(m["family"])
    assigned, unassigned = [], []
    for i, r in enumerate(rows):
        j = dict(contig=r["contig"], bp=r["bp"], side="cigar_insertion", families=sorted(fams.get(i, ())),
                 molecule=molecule(r["read"]), read=r["read"], length=r["length"])
        (assigned if j["families"] else unassigned).append(j)
    return assigned, unassigned


def cluster_sites(junctions, keep):
    rows = sorted((j for j in junctions if j["molecule"] in keep), key=lambda j: (j["contig"], j["bp"]))
    sites = []
    for j in rows:
        if sites and sites[-1]["contig"] == j["contig"] and j["bp"] - sites[-1]["end"] <= CLUSTER_BP:
            site = sites[-1]
        else:
            site = dict(contig=j["contig"], start=j["bp"], end=j["bp"], molecules=set(), sides=set(), families=set(),
                        lengths=[])
            sites.append(site)
        site["end"] = j["bp"]
        site["molecules"].add(j["molecule"])
        site["sides"].add(j["side"])
        site["families"].update(j["families"])
        if "length" in j:
            site["lengths"].append(j["length"])
    return sites


def two_sided(site):
    return "cigar_insertion" in site["sides"] or {"anchor_left", "anchor_right"} <= site["sides"]


def write_sites(path, sites):
    with open(path, "w") as out:
        out.write("contig\tstart0\tend0\tmolecules\ttwo_sided\tevidence\tmedian_insert_len\tfamilies\tcalled\n")
        for s in sites:
            lens = sorted(s["lengths"])
            out.write("\t".join(map(str, [s["contig"], s["start"], s["end"], len(s["molecules"]), two_sided(s),
                                          ",".join(sorted(s["sides"])), lens[len(lens) // 2] if lens else "",
                                          ",".join(sorted(s["families"])) or "-",
                                          len(s["molecules"]) >= MIN_SITE_MOLECULES])) + "\n")


def keep_fraction(mol, seed, fraction):
    h = int(hashlib.blake2b(("%d:%s" % (seed, mol)).encode(), digest_size=8).hexdigest(), 16)
    return h / 2**64 < fraction


def classify_main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    members = load_members(args.run)
    index = MemberIndex(members)
    contig_len = fasta_lengths(args.reference)
    mol_len, by_read = collect(args.paf, index)
    genome = sum(contig_len.values())
    depth = sum(mol_len.values()) / genome if genome else 0
    ev = member_evidence(by_read, index, contig_len)
    split_junctions, known_splits, palindromic = nonreference_junctions(by_read, index)
    ins_assigned, ins_unassigned = insertion_junctions(args.insertions, args.ins_paf, index)
    # A CIGAR insertion at a reference copy's boundary is not a non-reference site.
    ins_assigned = [j for j in ins_assigned if not any(
        m["family"] in j["families"] and min(abs(j["bp"] - m["s"]), abs(j["bp"] - m["e"])) <= JUNCTION_SLOP_BP
        for m in index.overlapping(j["contig"], j["bp"], j["bp"] + 1, pad=JUNCTION_SLOP_BP))]
    junctions = split_junctions + ins_assigned
    all_mols = set(mol_len)

    classes = {}
    with open(args.out / "member_read_support.tsv", "w") as out:
        out.write("member\tfamily\tmethod\ttnp_supported\tcontig\tstart0\tend0\tlength\tboundary_source\t"
                  "full\tleft\tright\tempty\tambiguous\tinternal\tclass\n")
        for m in members:
            mols = ev.get(m["member"], {})
            n = {f: sum(1 for x in mols.values() if f in x) for f in ("full", "left", "right", "empty", "ambiguous", "internal")}
            cls = classify(m, mols, contig_len)
            classes[m["member"]] = cls
            out.write("\t".join(map(str, [m["member"], m["family"], m["method"], m["tnp_supported"], m["contig"], m["s"], m["e"],
                                          m["e"] - m["s"], m["boundary_source"], n["full"], n["left"], n["right"], n["empty"],
                                          n["ambiguous"], n["internal"], cls])) + "\n")

    sites = cluster_sites(junctions, all_mols)
    write_sites(args.out / "nonreference_family_sites.tsv", sites)
    unassigned_sites = cluster_sites(ins_unassigned, all_mols)
    write_sites(args.out / "unassigned_insertion_sites.tsv", unassigned_sites)

    def called(sites_):
        return [x for x in sites_ if len(x["molecules"]) >= MIN_SITE_MOLECULES]

    curve = []
    targets = [d for d in DEPTHS if d < depth] + [round(depth, 1)]
    for target in targets:
        frac = target / depth if target in DEPTHS and depth else 1.0
        for seed in SEEDS:
            keep = {m for m in all_mols if keep_fraction(m, seed, frac)}
            counts = dict.fromkeys(CLASS_ORDER, 0)
            for m in members:
                mols = {k: v for k, v in ev.get(m["member"], {}).items() if k in keep}
                counts[classify(m, mols, contig_len)] += 1
            fam = called(cluster_sites(junctions, keep))
            una = called(cluster_sites(ins_unassigned, keep))
            curve.append(dict(depth=target, seed=seed, molecules=len(keep), **counts,
                              nonreference_family_sites=len(fam),
                              nonreference_family_two_sided=sum(two_sided(x) for x in fam),
                              unassigned_insertion_sites=len(una)))
            if frac >= 1.0:
                break
    with open(args.out / "depth_curve.tsv", "w") as out:
        cols = list(curve[0]) if curve else []
        out.write("\t".join(cols) + "\n")
        for row in curve:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")

    class_counts = dict.fromkeys(CLASS_ORDER, 0)
    for c in classes.values():
        class_counts[c] += 1
    fam_called, una_called = called(sites), called(unassigned_sites)
    summary = dict(
        sample=args.sample,
        molecules_mapped=len(all_mols),
        genome_bp=genome,
        molecule_depth=round(depth, 1),
        members=len(members),
        member_classes=class_counts,
        member_classes_tnp_supported={c: sum(1 for m in members if m["tnp_supported"] and classes[m["member"]] == c)
                                      for c in CLASS_ORDER},
        reads_touching_candidates=len(by_read),
        palindromic_read_pieces=palindromic,
        known_boundary_split_alignments=known_splits,
        split_read_nonreference_junctions=len(split_junctions),
        cigar_insertions_assigned_to_family=len(ins_assigned),
        cigar_insertions_unassigned=len(ins_unassigned),
        nonreference_family_sites_called=len(fam_called),
        nonreference_family_sites_two_sided=sum(two_sided(x) for x in fam_called),
        unassigned_insertion_sites_called=len(una_called),
        thresholds=dict(ANCHOR_BP=ANCHOR_BP, MIN_MAPQ=MIN_MAPQ, BIG_INDEL_FRACTION=BIG_INDEL_FRACTION,
                        EMPTY_GAP_BP=EMPTY_GAP_BP, JUNCTION_SLOP_BP=JUNCTION_SLOP_BP, CLUSTER_BP=CLUSTER_BP,
                        ADJACENT_READ_GAP=ADJACENT_READ_GAP, MIN_ANCHOR_LEN=MIN_ANCHOR_LEN,
                        MIN_SITE_MOLECULES=MIN_SITE_MOLECULES, MIN_INSERTION_BP=MIN_INSERTION_BP,
                        MIN_INSERT_COVERED=MIN_INSERT_COVERED, MIN_MEMBER_CONTAINED=MIN_MEMBER_CONTAINED, DEPTHS=DEPTHS, SEEDS=SEEDS),
    )
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("insertions")
    p.add_argument("--paf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("extract")
    p.add_argument("--insertions", type=Path, required=True)
    p.add_argument("--fastq", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("classify")
    p.add_argument("--run", type=Path, required=True, help="multisite-insertion run directory")
    p.add_argument("--reference", type=Path, required=True, help="assembly FASTA used for mapping")
    p.add_argument("--paf", type=Path, required=True)
    p.add_argument("--insertions", type=Path, required=True)
    p.add_argument("--ins-paf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--sample", required=True)
    args = parser.parse_args()
    if args.cmd == "insertions":
        scan_insertions(args.paf, args.out)
    elif args.cmd == "extract":
        extract_insertions(args.insertions, args.fastq, args.out)
    else:
        classify_main(args)


if __name__ == "__main__":
    main()
