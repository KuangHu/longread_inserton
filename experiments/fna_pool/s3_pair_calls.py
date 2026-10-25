"""S3 prototype: pairwise assembly gap calls (U-X-D in one genome vs U-D in the other), locus merging, U/core/D export.

call    For one pair, maps each genome onto the other (minimap2 -x asm20 -c --eqx, secondaries kept) and records
        every difference of >= MIN_CORE bp with two flanks: (a) I/D runs inside one alignment, (b) a gap between two
        collinear alignment blocks of the same contig pair and strand. Each raw call keeps both genomes' coordinates,
        mapping direction, alignment strand, flank identity/cover, contig-end, tandem-like, boundary shift window and
        other-side gap (replacement). Probe step: the flanks and core of every call are mapped as short queries to both
        genomes (a whole-genome alignment hides short secondary placements): a flank with > 1 placement in either genome
        is not unique; core copies in each genome are counted. Alignment ends where both genomes continue unaligned
        (e.g. next to a contig end) are written to one_flank.tsv as unresolved one-flank records, not loci.
        Presence/absence is recorded as "X present in carrier, absent at the empty site"; no insertion/deletion
        direction is inferred.
merge   Merges raw calls of all pairs into per-genome loci (same carrier contig; reciprocal overlap >= 0.8 or both
        boundaries within the shift window + BOUNDARY_TOL). Adjacent but non-overlapping cores stay separate.
        Writes loci.tsv, calls_merged.tsv and U/core/D FASTA.
plot    SVG per locus: carrier window with GFF3 CDS arrows, U/core/D marked, and the partner's empty site window.

Statuses: two_flank_anchored (both flanks >= FLANK_MIN_ID identity and >= FLANK_MIN_COVER aligned, no competing
placement of either flank, not within EDGE bp of a contig end, other-side gap small), replacement, or unresolved
with reasons. Thresholds are prototype settings for the first pool, not calibrated values.
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

VERSION = "s3_pair_calls_v0.1"
MM = ["-x", "asm20", "-c", "--eqx", "--secondary=yes", "-N", "20", "-p", "0.5"]
MIN_CORE = 50
MAX_CORE = 100000
FLANK = 200
FLANK_MIN_ID = 0.90
FLANK_MIN_COVER = 0.90
MERGE_MATCH = 5  # I/D runs separated by <= this many aligned bases form one event
BLOCK_OVERLAP = 50  # adjacent blocks may overlap this much on either genome (microhomology / alignment ends)
SMALL_OTHER_GAP = 10  # other-side gap allowed for a clean presence/absence call (also 10% of core)
EDGE = 1000
COMPETE_COVER = 0.80
COMPETE_ID = 0.98  # competing alignment identity >= this fraction of the flank identity
TANDEM_K, TANDEM_CONTAIN = 13, 0.5
BOUNDARY_TOL = 20
UD = 1000
PROBE_MM = ["-x", "map-pb", "-c", "--eqx", "--secondary=yes", "-N", "100", "-p", "0.1"]  # asm*/map-hifi presets suppress equal-scoring secondaries of short queries; map-pb reports them
FLANK_HIT_ID, FLANK_HIT_COVER = 0.90, 0.80  # a flank placement; > 1 in either genome = not unique
CORE_HIT_ID, CORE_HIT_COVER = 0.80, 0.50  # a homologous copy of the core
PROBE_CORE_MAX = 50000
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(str(path), "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def parse_paf(path):
    out = []
    for line in open(str(path)):
        f = line.rstrip("\n").split("\t")
        tags = dict((x[:2], x[5:]) for x in f[12:])
        ops = [(int(n), o) for n, o in CIGAR.findall(tags.get("cg", ""))]
        eq = sum(n for n, o in ops if o == "=")
        tot = sum(n for n, o in ops if o in "=XID")
        out.append(dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], tlen=int(f[6]),
                        ts=int(f[7]), te=int(f[8]), mapq=int(f[11]), tp=tags.get("tp", "P"), ops=ops,
                        identity=eq / tot if tot else 0.0, AS=int(tags.get("AS", 0))))
    return out


def walk(a):
    """Yield (op, n, rpos, qpos_al) along the alignment; qpos_al counts query bases in alignment orientation."""
    r, q = a["ts"], 0
    for n, o in a["ops"]:
        yield o, n, r, q
        if o in "=XMD":
            r += n
        if o in "=XMI":
            q += n


def q_fwd(a, q_al):
    """Alignment-oriented query offset -> forward query coordinate (position before which the offset sits)."""
    return a["qs"] + q_al if a["strand"] == "+" else a["qe"] - q_al


def ref_to_q(a, rpos):
    for o, n, r, q in walk(a):
        if o in "=XMD" and r <= rpos <= r + n:
            return q_fwd(a, q + (rpos - r if o != "D" else 0))
    return None


def window_stats(a, rs, re_):
    """Identity and aligned cover of alignment a over reference window [rs, re_)."""
    eq = x = ins = dl = 0
    for o, n, r, q in walk(a):
        if o in "=XM":
            ov = max(0, min(r + n, re_) - max(r, rs))
            if o == "=":
                eq += ov
            else:
                x += ov
        elif o == "D":
            dl += max(0, min(r + n, re_) - max(r, rs))
        elif o == "I" and rs < r < re_:
            ins += n
    tot = eq + x + ins + dl
    return (eq / tot if tot else 0.0), (eq + x) / max(1, re_ - rs)


def competing(a, alns, qs, qe, flank_id):
    """Another alignment placing query window [qs, qe) elsewhere with comparable identity."""
    if qs is None or qe is None or qe <= qs:
        return []
    out = []
    for b in alns:
        if b is a or b["q"] != a["q"]:
            continue
        cov = max(0, min(b["qe"], qe) - max(b["qs"], qs)) / (qe - qs)
        same_place = b["t"] == a["t"] and b["ts"] < a["te"] and b["te"] > a["ts"]
        if cov >= COMPETE_COVER and not same_place and b["identity"] >= COMPETE_ID * flank_id:
            out.append("%s:%d-%d(%s,%.3f)" % (b["t"], b["ts"], b["te"], b["tp"], b["identity"]))
    return out


def shift_window(seq, s, e):
    """How far an indel of seq[s:e] can slide left/right with the same result (repeat-ambiguous boundary)."""
    left = 0
    while s - left - 1 >= 0 and seq[s - left - 1] == seq[e - left - 1]:
        left += 1
    right = 0
    while e + right < len(seq) and seq[s + right] == seq[e + right]:
        right += 1
    return left, right


def tandem_like(seq, s, e):
    L = e - s
    if L < TANDEM_K * 2:
        return False
    core = {seq[i:i + TANDEM_K] for i in range(s, e - TANDEM_K + 1)}
    for a, b in ((max(0, s - L), s), (e, min(len(seq), e + L))):
        side = {seq[i:i + TANDEM_K] for i in range(a, b - TANDEM_K + 1)}
        if core and len(core & side) / len(core) >= TANDEM_CONTAIN:
            return True
    return False


def events_in_alignment(a):
    """I/D runs >= MIN_CORE inside one alignment: (ref_s, ref_e, q_al_s, q_al_e)."""
    ev, cur, gap_since = [], None, 0
    for o, n, r, q in walk(a):
        if o in "ID":
            if cur is None:
                cur = [r, r, q, q]
            if o == "D":
                cur[1] = r + n
            else:
                cur[3] = q + n
            gap_since = 0
        elif o in "=XM":
            if cur is not None and n > MERGE_MATCH:
                ev.append(tuple(cur))
                cur = None
            elif cur is not None:
                cur[1], cur[3] = r + n, q + n
    if cur is not None:
        ev.append(tuple(cur))
    return [e for e in ev if max(e[1] - e[0], e[3] - e[2]) >= MIN_CORE]


def block_events(alns):
    """Gaps between collinear primary blocks of the same contig pair and strand."""
    groups = defaultdict(list)
    for a in alns:
        if a["tp"] == "P":
            groups[(a["t"], a["q"], a["strand"])].append(a)
    out = []
    for (t, q, strand), bl in groups.items():
        bl.sort(key=lambda a: a["ts"])
        for b1, b2 in zip(bl, bl[1:]):
            rgap = b2["ts"] - b1["te"]
            qgap = (b2["qs"] - b1["qe"]) if strand == "+" else (b1["qs"] - b2["qe"])
            if min(rgap, qgap) < -BLOCK_OVERLAP or max(rgap, qgap) < MIN_CORE or max(rgap, qgap) > MAX_CORE:
                continue
            out.append((b1, b2, rgap, qgap))
    return out


def call_pair(pair_id, ref_name, ref_seqs, qry_name, qry_seqs, paf, fastas=None, workdir=None, threads=2):
    """fastas: name -> FASTA path for the probe step (flank uniqueness and core copies in both genomes)."""
    alns = parse_paf(paf)
    calls, one_flank = [], []

    def record(source, a_left, a_right, rs, re_, qs, qe):
        """rs..re_ ref interval of the difference; qs..qe forward query interval (both may be empty)."""
        rlen, qlen = re_ - rs, qe - qs
        if max(rlen, qlen) < MIN_CORE or max(rlen, qlen) > MAX_CORE:
            return
        x_in_ref = rlen >= qlen
        l_id, l_cov = window_stats(a_left, rs - FLANK, rs)
        r_id, r_cov = window_stats(a_right, re_, re_ + FLANK)
        lq = sorted([ref_to_q(a_left, rs - FLANK), ref_to_q(a_left, rs)], key=lambda v: (v is None, v))
        rq = sorted([ref_to_q(a_right, re_), ref_to_q(a_right, re_ + FLANK)], key=lambda v: (v is None, v))
        comp = competing(a_left, alns, lq[0], lq[1], l_id) + competing(a_right, alns, rq[0], rq[1], r_id)  # within-mapping
        rseq, qseq = ref_seqs[a_left["t"]], qry_seqs[a_left["q"]]
        if x_in_ref:
            carrier, ccontig, cs, ce, cseq = ref_name, a_left["t"], rs, re_, rseq
            empty, econtig, epos, eseq = qry_name, a_left["q"], qs, qseq
            other_gap = qlen
        else:
            carrier, ccontig, cs, ce, cseq = qry_name, a_left["q"], qs, qe, qseq
            empty, econtig, epos, eseq = ref_name, a_left["t"], rs, rseq
            other_gap = rlen
        shl, shr = shift_window(cseq, cs, ce)
        near_edge = cs < EDGE or len(cseq) - ce < EDGE or epos < EDGE or len(eseq) - epos < EDGE
        reasons = []
        if min(l_id, r_id) < FLANK_MIN_ID or min(l_cov, r_cov) < FLANK_MIN_COVER:
            reasons.append("flank_identity_or_cover")
        if comp:
            reasons.append("flank_competing")
        if near_edge:
            reasons.append("near_contig_end")
        replacement = other_gap > max(SMALL_OTHER_GAP, 0.1 * (ce - cs))
        status = "replacement" if replacement and not reasons else ("two_flank_anchored" if not reasons else "unresolved")
        calls.append(dict(
            pair_id=pair_id, mapping="%s_on_%s" % (qry_name, ref_name), source=source, aln_strand=a_left["strand"],
            carrier=carrier, carrier_contig=ccontig, core_start=cs, core_end=ce, core_len=ce - cs,
            empty=empty, empty_contig=econtig, empty_pos=epos, other_side_gap=other_gap,
            flank_left_id=round(l_id, 4), flank_left_cover=round(l_cov, 4), flank_right_id=round(r_id, 4),
            flank_right_cover=round(r_cov, 4), mapq_left=a_left["mapq"], mapq_right=a_right["mapq"],
            flank_competing=";".join(comp), near_contig_end=near_edge, tandem_like=tandem_like(cseq, cs, ce),
            shift_left=shl, shift_right=shr, status=status, reasons=";".join(reasons), caller=VERSION))

    for a in alns:
        if a["tp"] != "P":
            continue
        for rs, re_, qa, qb in events_in_alignment(a):
            q1, q2 = sorted((q_fwd(a, qa), q_fwd(a, qb)))
            record("cigar", a, a, rs, re_, q1, q2)
    for b1, b2, rgap, qgap in block_events(alns):
        rs, re_ = b1["te"], max(b1["te"], b2["ts"])
        if b1["strand"] == "+":
            q1, q2 = b1["qe"], max(b1["qe"], b2["qs"])
        else:
            q1, q2 = b2["qe"], max(b2["qe"], b1["qs"])
        record("split_blocks", b1, b2, rs, re_, q1, q2)
    one_flank = alignment_end_events(pair_id, ref_name, ref_seqs, qry_name, qry_seqs, alns)
    if fastas and calls:
        probe(calls, {ref_name: ref_seqs, qry_name: qry_seqs}, fastas, workdir, threads)
    return calls, one_flank


def alignment_end_events(pair_id, ref_name, ref_seqs, qry_name, qry_seqs, alns):
    """Primary alignment ends where both genomes continue with >= MIN_CORE unaligned bases and no collinear block
    of the same contig pair follows: one-flank context only (e.g. a difference next to a contig end). Never anchored."""
    prim = [a for a in alns if a["tp"] == "P"]
    out = []
    for a in prim:
        for side in ("left", "right"):
            if side == "left":
                r_rest = a["ts"]
                q_rest = a["qs"] if a["strand"] == "+" else a["qlen"] - a["qe"]
                covered = any(b is not a and b["t"] == a["t"] and b["q"] == a["q"] and b["strand"] == a["strand"]
                              and b["te"] <= a["ts"] + BLOCK_OVERLAP and a["ts"] - b["te"] <= MAX_CORE for b in prim)
            else:
                r_rest = a["tlen"] - a["te"]
                q_rest = a["qlen"] - a["qe"] if a["strand"] == "+" else a["qs"]
                covered = any(b is not a and b["t"] == a["t"] and b["q"] == a["q"] and b["strand"] == a["strand"]
                              and b["ts"] >= a["te"] - BLOCK_OVERLAP and b["ts"] - a["te"] <= MAX_CORE for b in prim)
            if covered or min(r_rest, q_rest) <= 0 or max(r_rest, q_rest) < MIN_CORE:
                continue
            rpos = a["ts"] if side == "left" else a["te"]
            qpos = ref_to_q(a, rpos)
            out.append(dict(pair_id=pair_id, mapping="%s_on_%s" % (qry_name, ref_name), ref=ref_name, ref_contig=a["t"],
                            ref_pos=rpos, ref_unaligned_beyond=r_rest, qry=qry_name, qry_contig=a["q"], qry_pos=qpos,
                            qry_unaligned_beyond=q_rest, aln_strand=a["strand"], side=side,
                            ref_contig_end_within=min(rpos, a["tlen"] - rpos), qry_contig_end_within=min(qpos, a["qlen"] - qpos),
                            near_contig_end=min(rpos, a["tlen"] - rpos) < EDGE or min(qpos, a["qlen"] - qpos) < EDGE,
                            status="unresolved", reasons="one_flank_only", caller=VERSION))
    return out


def probe(calls, seqs, fastas, workdir, threads):
    """Map each call's flanks and core (short queries, so every copy scores comparably) to both genomes."""
    workdir = Path(workdir)
    q = workdir / "probes.fa"
    with open(str(q), "w") as h:
        for i, c in enumerate(calls):
            cs = seqs[c["carrier"]][c["carrier_contig"]]
            s, e = c["core_start"], c["core_end"]
            for part, ps, pe in (("L", max(0, s - FLANK), s), ("R", e, min(len(cs), e + FLANK)), ("C", s, min(e, s + PROBE_CORE_MAX))):
                if pe - ps >= 30:
                    h.write(">%d|%s\n%s\n" % (i, part, cs[ps:pe]))
    hits = defaultdict(Counter)
    for name, fa in fastas.items():
        paf = workdir / ("probes_on_%s.paf" % name)
        with open(str(paf), "w") as h:
            subprocess.run(["minimap2"] + PROBE_MM + ["-t", str(threads), str(fa), str(q)], stdout=h,
                           stderr=subprocess.DEVNULL, check=True)
        for b in parse_paf(paf):
            i, part = b["q"].split("|")
            cover = (b["qe"] - b["qs"]) / b["qlen"]
            if part in "LR" and b["identity"] >= FLANK_HIT_ID and cover >= FLANK_HIT_COVER:
                hits[(int(i), name)][part] += 1
            if part == "C" and b["identity"] >= CORE_HIT_ID and cover >= CORE_HIT_COVER:
                hits[(int(i), name)]["C"] += 1
    for i, c in enumerate(calls):
        hc, he = hits[(i, c["carrier"])], hits[(i, c["empty"])]
        c.update(flank_hits_carrier="L%d/R%d" % (hc["L"], hc["R"]), flank_hits_empty="L%d/R%d" % (he["L"], he["R"]),
                 core_copies_carrier=hc["C"], core_copies_empty=he["C"])
        multi = [x for x, n in (("L_carrier", hc["L"]), ("R_carrier", hc["R"]), ("L_empty", he["L"]), ("R_empty", he["R"])) if n > 1]
        if multi:
            reasons = [r for r in c["reasons"].split(";") if r] + ["flank_not_unique:" + ",".join(multi)]
            c["reasons"] = ";".join(reasons)
            c["status"] = "unresolved"


def run_minimap2(ref, qry, out, threads):
    with open(str(out), "w") as h:
        subprocess.run(["minimap2"] + MM + ["-t", str(threads), str(ref), str(qry)], stdout=h, stderr=subprocess.DEVNULL, check=True)


CALL_COLS = None


def cmd_call(a):
    a.out.mkdir(parents=True, exist_ok=True)
    seqs = {a.name_a: read_fasta(a.fasta_a), a.name_b: read_fasta(a.fasta_b)}
    calls, one = [], []
    fastas = {a.name_a: a.fasta_a, a.name_b: a.fasta_b}
    for ref, qry, fr, fq in ((a.name_a, a.name_b, a.fasta_a, a.fasta_b), (a.name_b, a.name_a, a.fasta_b, a.fasta_a)):
        paf = a.out / ("%s_on_%s.paf" % (qry, ref))
        run_minimap2(fr, fq, paf, a.threads)
        wd = a.out / ("probe_%s_on_%s" % (qry, ref))
        wd.mkdir(exist_ok=True)
        c, o = call_pair(a.pair_id, ref, seqs[ref], qry, seqs[qry], paf, fastas, wd, a.threads)
        calls, one = calls + c, one + o
    write(a.out / "calls.tsv", calls)
    write(a.out / "one_flank.tsv", one)
    ver = subprocess.run(["minimap2", "--version"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    (a.out / "call_run.json").write_text(json.dumps(dict(pair_id=a.pair_id, caller=VERSION, minimap2=ver, minimap2_args=MM,
                                                         calls=len(calls), one_flank=len(one),
                                                         status=dict(Counter(c["status"] for c in calls))), indent=2) + "\n")
    print(a.pair_id, len(calls), dict(Counter(c["status"] for c in calls)))


def write(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else ["pair_id"])
    with open(str(path), "w") as h:
        w = csv.DictWriter(h, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def merge_calls(calls):
    """Per carrier contig: same locus if reciprocal overlap >= 0.8 or both boundaries within shift window + tol."""
    by = defaultdict(list)
    for c in calls:
        by[(c["carrier"], c["carrier_contig"])].append(c)
    loci = []
    for (carrier, contig), cs in sorted(by.items()):
        cs.sort(key=lambda c: (int(c["core_start"]), int(c["core_end"])))
        clusters = []
        for c in cs:
            s, e = int(c["core_start"]), int(c["core_end"])
            tol = BOUNDARY_TOL + max(int(c["shift_left"]), int(c["shift_right"]))
            for cl in clusters:
                rs, re_ = cl["s"], cl["e"]
                ov = max(0, min(e, re_) - max(s, rs))
                if (ov >= 0.8 * (e - s) and ov >= 0.8 * (re_ - rs)) or (abs(s - rs) <= tol and abs(e - re_) <= tol):
                    cl["calls"].append(c)
                    break
            else:
                clusters.append(dict(s=s, e=e, calls=[c]))
        for cl in clusters:
            good = [c for c in cl["calls"] if c["status"] == "two_flank_anchored"] or cl["calls"]
            (s, e), _ = Counter((int(c["core_start"]), int(c["core_end"])) for c in good).most_common(1)[0]
            lid = "L_" + hashlib.sha256(("%s|%s|%d|%d" % (carrier, contig, s, e)).encode()).hexdigest()[:12]
            st = Counter(c["status"] for c in cl["calls"])
            loci.append(dict(
                locus_id=lid, carrier=carrier, contig=contig, core_start=s, core_end=e, core_len=e - s,
                status="two_flank_anchored" if st["two_flank_anchored"] else ("replacement" if st["replacement"] else "unresolved"),
                n_calls=len(cl["calls"]), n_pairs=len({c["pair_id"] for c in cl["calls"]}),
                empty_partners=";".join(sorted({c["empty"] for c in cl["calls"]})),
                mappings=len({(c["pair_id"], c["mapping"]) for c in cl["calls"]}), call_status=json.dumps(dict(st)),
                start_spread=max(int(c["core_start"]) for c in cl["calls"]) - min(int(c["core_start"]) for c in cl["calls"]),
                end_spread=max(int(c["core_end"]) for c in cl["calls"]) - min(int(c["core_end"]) for c in cl["calls"]),
                boundary_consistency=boundary_consistency(cl["calls"]),
                tandem_like=any(c["tandem_like"] in (True, "True") for c in cl["calls"]),
                core_copies_carrier=max(int(c.get("core_copies_carrier") or 0) for c in cl["calls"]),
                core_copies_empty=max(int(c.get("core_copies_empty") or 0) for c in cl["calls"]),
                near_contig_end=any(c["near_contig_end"] in (True, "True") for c in cl["calls"]),
                reasons=";".join(sorted({r for c in cl["calls"] for r in c["reasons"].split(";") if r}))))
            for c in cl["calls"]:
                c["locus_id"] = lid
    return loci


def boundary_consistency(calls):
    """exact: all calls give the same core; within_shift: differences fit the indel shift window + BOUNDARY_TOL;
    partner_dependent: calls against different partners give different boundaries (raw calls are kept)."""
    ss = {int(c["core_start"]) for c in calls}
    ee = {int(c["core_end"]) for c in calls}
    if len(ss) == 1 and len(ee) == 1:
        return "exact"
    tol = BOUNDARY_TOL + max(max(int(c["shift_left"]), int(c["shift_right"])) for c in calls)
    if max(ss) - min(ss) <= tol and max(ee) - min(ee) <= tol:
        return "within_shift"
    return "partner_dependent"


def cmd_merge(a):
    calls = []
    for p in a.calls:
        calls += list(csv.DictReader(open(str(p)), delimiter="\t"))
    loci = merge_calls(calls)
    a.out.mkdir(parents=True, exist_ok=True)
    write(a.out / "calls_merged.tsv", calls)
    write(a.out / "loci.tsv", loci)
    if a.fasta:
        fa = dict(x.split("=", 1) for x in a.fasta)
        seqs = {k: read_fasta(v) for k, v in fa.items()}
        with open(str(a.out / "ucd.fasta"), "w") as h:
            for l in loci:
                seq = seqs[l["carrier"]][l["contig"]]
                s, e = l["core_start"], l["core_end"]
                for part, ps, pe in (("U", max(0, s - UD), s), ("core", s, e), ("D", e, min(len(seq), e + UD))):
                    h.write(">%s|%s|%s:%d-%d|%s\n%s\n" % (l["locus_id"], l["carrier"], l["contig"], ps, pe, part, seq[ps:pe]))
    print(len(calls), "calls ->", len(loci), "loci", dict(Counter(l["status"] for l in loci)))


def read_gff_cds(path):
    cds = defaultdict(list)
    if not path:
        return cds
    for line in open(str(path)):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("CDS", "rRNA", "tRNA"):
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        cds[f[0]].append((int(f[3]) - 1, int(f[4]), f[6], f[2], at.get("product", "")[:40], at.get("pseudo") == "true"))
    return cds


def svg_track(y, label, contig_len, ws, we, feats, marks, width=1000):
    scale = width / max(1, we - ws)
    x = lambda p: 60 + (p - ws) * scale
    out = ['<text x="4" y="%d" font-size="11">%s</text>' % (y - 14, label),
           '<line x1="60" y1="%d" x2="%d" y2="%d" stroke="#888"/>' % (y, 60 + width, y)]
    for s, e, color, name in marks:
        out.append('<rect x="%.1f" y="%d" width="%.1f" height="10" fill="%s" opacity="0.45"/><text x="%.1f" y="%d" font-size="9">%s</text>'
                   % (x(max(s, ws)), y - 5, max(1, x(min(e, we)) - x(max(s, ws))), color, x(max(s, ws)), y + 18, name))
    seen = set()
    for s, e, strand, typ, prod, pseudo in feats:
        if e < ws or s > we:
            continue
        x1, x2 = x(max(s, ws)), x(min(e, we))
        yy = y - 16 if strand == "+" else y + 4
        col = {"CDS": "#3b6fb6", "rRNA": "#b63b3b", "tRNA": "#b68a3b"}[typ] if not pseudo else "#999"
        if strand == "+":
            base = max(x1, x2 - 6)
            pts = [(x1, yy), (base, yy), (x2, yy + 6), (base, yy + 12), (x1, yy + 12)]
        else:
            base = min(x2, x1 + 6)
            pts = [(x1, yy + 6), (base, yy), (x2, yy), (x2, yy + 12), (base, yy + 12)]
        out.append('<polygon points="%s" fill="%s"><title>%s %d-%d %s</title></polygon>'
                   % (" ".join("%.1f,%d" % p for p in pts), col, prod, s, e, strand))
        key = (round(x1 / 40), strand)
        if prod and key not in seen and x2 - x1 > 30:
            seen.add(key)
            out.append('<text x="%.1f" y="%d" font-size="8" fill="#333">%s</text>' % (x1, yy - 2 if strand == "+" else yy + 22, prod[:24]))
    out.append('<text x="60" y="%d" font-size="9" fill="#666">%s:%d-%d</text>' % (y + 34, label.split()[0], ws, we))
    return out


def cmd_plot(a):
    loci = list(csv.DictReader(open(str(a.loci)), delimiter="\t"))
    calls = list(csv.DictReader(open(str(a.calls)), delimiter="\t"))
    fa = dict(x.split("=", 1) for x in a.fasta)
    gff = dict(x.split("=", 1) for x in (a.gff or []))
    lens = {k: {c: len(s) for c, s in read_fasta(v).items()} for k, v in fa.items()}
    feats = {k: read_gff_cds(v) for k, v in gff.items()}
    a.out.mkdir(parents=True, exist_ok=True)
    want = set(a.loci_ids or [])
    for l in loci:
        if want and l["locus_id"] not in want:
            continue
        s, e = int(l["core_start"]), int(l["core_end"])
        ws, we = max(0, s - UD - a.pad), min(lens[l["carrier"]][l["contig"]], e + UD + a.pad)
        c = [c for c in calls if c.get("locus_id") == l["locus_id"]][0]
        ep = int(c["empty_pos"])
        es, ee = max(0, ep - UD - a.pad), min(lens[c["empty"]][c["empty_contig"]], ep + UD + a.pad + (e - s) // 4)
        body = ['<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="250" font-family="sans-serif">',
                '<rect width="100%" height="100%" fill="white"/>',
                '<text x="4" y="14" font-size="12">%s  %s  core %d bp  status %s  pairs %s  tandem %s  core copies carrier/empty %s</text>'
                % (l["locus_id"], l["carrier"], e - s, l["status"], l["n_pairs"], l["tandem_like"],
                   "%s/%s" % (l["core_copies_carrier"], l["core_copies_empty"]))]
        body += svg_track(80, "%s %s (carrier)" % (l["carrier"], l["contig"]), 0, ws, we, feats.get(l["carrier"], {}).get(l["contig"], []),
                          [(max(0, s - UD), s, "#2a9d8f", "U"), (s, e, "#e76f51", "core"), (e, e + UD, "#2a9d8f", "D")])
        body += svg_track(190, "%s %s (empty site, %s)" % (c["empty"], c["empty_contig"], c["aln_strand"]), 0, es, ee,
                          feats.get(c["empty"], {}).get(c["empty_contig"], []), [(ep - 1, ep + 1, "#e76f51", "site")])
        body.append("</svg>")
        (a.out / (l["locus_id"] + ".svg")).write_text("\n".join(body) + "\n")
    print("plots in", a.out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode")
    c = sub.add_parser("call")
    c.add_argument("--pair-id", required=True)
    c.add_argument("--name-a", required=True)
    c.add_argument("--fasta-a", type=Path, required=True)
    c.add_argument("--name-b", required=True)
    c.add_argument("--fasta-b", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--threads", type=int, default=4)
    m = sub.add_parser("merge")
    m.add_argument("--calls", type=Path, nargs="+", required=True)
    m.add_argument("--fasta", nargs="*", help="name=path for U/core/D export")
    m.add_argument("--out", type=Path, required=True)
    g = sub.add_parser("plot")
    g.add_argument("--loci", type=Path, required=True)
    g.add_argument("--calls", type=Path, required=True)
    g.add_argument("--fasta", nargs="+", required=True)
    g.add_argument("--gff", nargs="*")
    g.add_argument("--loci-ids", nargs="*")
    g.add_argument("--pad", type=int, default=2000)
    g.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    dict(call=cmd_call, merge=cmd_merge, plot=cmd_plot)[a.mode](a)


if __name__ == "__main__":
    main()
