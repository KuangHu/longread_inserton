"""S3 long-read allele check: do a sample's own long reads support the allele its assembly shows at a call?

Unit = one test: a distinct raw call (carrier core coordinates + partner empty site) of a selected locus. Each
sample is tested against two local alleles built from its own assembly, so flank divergence between strains cannot
favour one allele:
  carrier (role C): present = U + core + D, absent = U + other + D   (other = partner's other-side segment, oriented)
  partner (role E): absent = U' + other + D', present = U' + core + D' (core oriented to the partner's strand)
with FLANK bp of U/D. A read record supports an allele only when all of these hold:
  span      its alignment to that allele covers both junctions plus ANCHOR bp (plus the boundary shift window);
  CIGAR     inside that interval no merged indel event (runs separated by <= MERGE_MATCH aligned bases) and no net
            indel excess over the read's own net indel rate outside the interval reaches half the length difference
            between the two alleles (at least MIN_INDEL): a read of the other allele forced onto this one shows that
            difference as a deletion or insertion, so an alignment that deletes the core cannot count as present;
  allele    if it passes on both alleles, the higher alignment score wins by >= MIN_DELTA_AS, else ambiguous;
  locus     the read is also mapped to the sample's whole assembly. The two read positions that the supporting
            local alignment puts at the outer ends of the tested interval (U side and D side anchors) must each be
            placed by a genome alignment of the same strand within LOCUS_TOL bp of where that allele position sits in
            the assembly, and no genome alignment placing the same read position elsewhere may score within
            MIN_DELTA_AS of it (off_locus / locus_ambiguous otherwise). Repeat copies thus compete genome-wide on the
            read segment that actually supports the allele; touching the site window is not enough. A read counts for
            several loci only when its own anchors sit at each of them (true neighbouring sites). The anchors are
            outside the core, so a read of the other allele is not penalised for not matching the assembly's core.
Read records are then grouped into molecules: PacBio ZMW (movie/zmw) or SRA spot name; ONT read ids, with duplex
records ("id1;id2") joined to their simplex parents. Records of one molecule that disagree make it ambiguous. Runs
whose names carry no molecule information (molecule_grouping_unavailable) are counted in read records and labelled
so (count_unit); there ">= MIN_SPAN" means read records, not independent molecules.

Length bias: a longer allele needs longer reads to span. Verdicts use length-matched counts: a molecule counts only
if its read length >= the longer of the two tested intervals (raw counts are also reported). This removes the hard
cut-off; reads close to that length still span the shorter allele more often.
Role verdicts (expected allele: C present, E absent), from length-matched molecules n = expected + other:
  insufficient (n < MIN_SPAN), consistent (other <= MAX_OTHER_FRACTION * n), contradicted (expected <= that
  fraction), mixed (both alleles). Unchanged from v0.1; fixed in advance and not tuned on the pilot.
Test verdicts: both_samples_consistent, carrier_only_consistent (carrier reads consistent, partner without reads or
insufficient), partner_only_consistent (the reverse; the carrier assembly is then not checked), mixed, contradicted,
not_tested. Locus evidence_level: E2 all tests both samples, E2s some tests both samples, E1cp carrier and partner
consistent only in different tests, E1c carrier only, E1p partner only, E0 not tested, X mixed or contradicted.
The reads are those deposited for the assembly's BioSample and were usually used to build it: this is an
assembly-consistency check (misassembly, collapsed repeats, within-sample allele mixtures), not independent
sequencing validation.

alleles   tests.tsv + alleles_<sample>.fa from loci/calls (loci with status two_flank_anchored; all their calls).
fetch     Download the usable runs of the tested samples (scripts/fetch_long_reads.py fetch logic, MD5 checked).
map       minimap2 (preset by read type) of one sample's runs against its allele FASTA.
genome    For one sample: subset the reads that hit its alleles, record their molecule keys, map them to the
          sample's whole assembly (secondaries kept) for the locus check.
summarize Per test and per locus counts and verdicts.
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

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from s3_pair_calls import CIGAR, DIFF_ORDER, MERGE_MATCH, difference_class, read_fasta, rc  # noqa: E402

VERSION = "s3_readcheck_v0.3"
FLANK = 2000
ANCHOR = 300
MIN_DELTA_AS = 20
MIN_SPAN = 3
MAX_OTHER_FRACTION = 0.1
MIN_INDEL = 10  # floor of the CIGAR limit when the two alleles differ little in length
LOCUS_TOL = 50  # bp between expected and observed genome position of an anchor (< ANCHOR)
PRESET = {"ont": "map-ont", "pacbio_hifi_ccs": "map-hifi", "pacbio_clr_subreads": "map-pb", "pacbio_unresolved": "map-pb"}
MM = ["-c", "--secondary=yes", "-N", "100", "-p", "0.3"]
MM_GENOME = ["-c", "--secondary=yes", "-N", "50", "-p", "0.5"]
ONT_ID = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-")


def read(path):
    return list(csv.DictReader(open(str(path)), delimiter="\t"))


def write(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else ["test_id"])
    with open(str(path), "w") as h:
        w = csv.DictWriter(h, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def build_alleles(c, cseq, eseq):
    """Return {(role, allele): (seq, j1, j2, shift_left, shift_right)}; j1/j2 are the junctions in that allele."""
    cs, ce = int(c["core_start"]), int(c["core_end"])
    ep, gap = int(c["empty_pos"]), int(c["other_side_gap"])
    minus = c["aln_strand"] == "-"
    core, other = cseq[cs:ce], eseq[ep:ep + gap]
    shl, shr = int(c["shift_left"]), int(c["shift_right"])
    out = {}
    u, d = cseq[max(0, cs - FLANK):cs], cseq[ce:ce + FLANK]
    other_c = rc(other) if minus else other
    out[("C", "P")] = (u + core + d, len(u), len(u) + len(core), shl, shr)
    out[("C", "A")] = (u + other_c + d, len(u), len(u) + len(other_c), shl, shr)
    u2, d2 = eseq[max(0, ep - FLANK):ep], eseq[ep + gap:ep + gap + FLANK]
    core_e = rc(core) if minus else core
    sl, sr = (shr, shl) if minus else (shl, shr)
    out[("E", "A")] = (u2 + other + d2, len(u2), len(u2) + len(other), sl, sr)
    out[("E", "P")] = (u2 + core_e + d2, len(u2), len(u2) + len(core_e), sl, sr)
    return out


def cmd_alleles(a):
    keep = {x.split("\t")[0] for x in open(str(a.loci_ids)).read().split("\n")[1:] if x} if a.loci_ids else None
    loci = {l["locus_id"]: l for l in read(a.loci) if keep is None or l["locus_id"] in keep}
    fa = dict(x.split("=", 1) for x in a.fasta)
    seqs = {}
    tests, alleles = {}, defaultdict(list)
    for c in read(a.calls):
        l = loci.get(c["locus_id"])
        if l is None:
            continue
        if l["status"] != "two_flank_anchored":
            continue
        key = "|".join(c[k] for k in ("carrier", "carrier_contig", "core_start", "core_end", "empty", "empty_contig",
                                      "empty_pos", "other_side_gap", "aln_strand"))
        tid = "T_" + hashlib.sha256(key.encode()).hexdigest()[:12]
        if tid in tests:
            tests[tid]["calls"] += 1
            tests[tid]["call_status"][c["status"]] += 1
            continue
        for s in (c["carrier"], c["empty"]):
            if s not in seqs:
                seqs[s] = read_fasta(fa[s])
        cseq, eseq = seqs[c["carrier"]][c["carrier_contig"]], seqs[c["empty"]][c["empty_contig"]]
        al = build_alleles(c, cseq, eseq)
        cs, ce, ep, gap = int(c["core_start"]), int(c["core_end"]), int(c["empty_pos"]), int(c["other_side_gap"])
        other = eseq[ep:ep + gap]
        dclass, norm_gap = difference_class(cseq[cs:ce], rc(other) if c["aln_strand"] == "-" else other)
        t = dict(test_id=tid, locus_id=c["locus_id"], pair_id=c["pair_id"], carrier=c["carrier"], empty=c["empty"],
                 core_len=ce - cs, other_side_gap=gap, aln_strand=c["aln_strand"], call_status_first=c["status"],
                 locus_support=l.get("anchored_support", ""), boundary_consistency=l["boundary_consistency"],
                 difference_class=dclass, other_norm_gap=norm_gap,
                 C_site="%s:%d-%d" % (c["carrier_contig"], cs, ce), E_site="%s:%d-%d" % (c["empty_contig"], ep, ep + gap),
                 calls=1, call_status=Counter([c["status"]]))
        for (role, allele), (seq, j1, j2, sl, sr) in al.items():
            sample = c["carrier"] if role == "C" else c["empty"]
            name = "%s|%s|%s" % (tid, role, allele)
            alleles[sample].append((name, seq))
            t.update({"%s%s_len" % (role, allele): len(seq), "%s%s_j" % (role, allele): "%d,%d,%d,%d" % (j1, j2, sl, sr)})
        tests[tid] = t
    a.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for t in tests.values():
        t["call_status"] = json.dumps(dict(t["call_status"]))
        rows.append(t)
    write(a.out / "tests.tsv", rows)
    for s, al in alleles.items():
        with open(str(a.out / ("alleles_%s.fa" % s)), "w") as h:
            for name, seq in al:
                h.write(">%s\n%s\n" % (name, seq))
    print(len(rows), "tests;", len(alleles), "samples;", Counter(t["call_status_first"] for t in rows))


def sample_runs(pool_dir, samples):
    runs = {r["run_accession"]: r for r in read(pool_dir / "runs.tsv")}
    out = defaultdict(list)
    for r in read(pool_dir / "assembly_run_links.tsv"):
        if r["biosample"] in samples and r["run_usable"] == "True" and r["run"] not in {x["run_accession"] for x in out[r["biosample"]]}:
            out[r["biosample"]].append(runs[r["run"]])
    return out


def cmd_fetch(a):
    import fetch_long_reads as flr  # noqa: E402  (curl download, MD5 check, verify.json; refuses to overwrite)
    from concurrent.futures import ThreadPoolExecutor
    tests = read(a.work / "tests.tsv")
    samples = {t["carrier"] for t in tests} | {t["empty"] for t in tests}
    todo = []
    for s, rs in sorted(sample_runs(a.pool_dir, samples).items()):
        for r in rs:
            if (a.reads / r["run_accession"] / "verify.json").exists():
                print(s, r["run_accession"], "already fetched")
            elif ";" in r["fastq_ftp"]:
                print(s, r["run_accession"], "skipped: more than one FASTQ file", file=sys.stderr)
            else:
                todo.append((s, r))

    def one(item):  # runs go to separate directories, so several can download at once
        s, r = item
        run = dict(run=r["run_accession"], sample_id=s, url="https://" + r["fastq_ftp"], fastq_bytes=int(r["fastq_bytes"]),
                   fastq_md5=r["fastq_md5"], platform=r["instrument_platform"], archive_base_count=int(r["base_count"]),
                   archive_read_count=int(r["read_count"]))
        try:
            return s, r["run_accession"], flr.fetch(run, a.reads)["status"]
        except Exception as e:  # reported, not fatal: map uses only verified runs
            return s, r["run_accession"], "error: %r" % e

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for s, run, st in ex.map(one, todo):
            print(s, run, st, flush=True)


def verified_runs(a):
    out = []
    for r in sample_runs(a.pool_dir, {a.sample})[a.sample]:
        v = a.reads / r["run_accession"] / "verify.json"
        if not v.exists():
            print(a.sample, r["run_accession"], "not fetched", file=sys.stderr)
            continue
        v = json.load(open(str(v)))
        if v["status"] != "verified":
            print(a.sample, r["run_accession"], "not verified:", v["checks"], file=sys.stderr)
            continue
        out.append((r, v))
    return out


def cmd_map(a):
    fa = a.work / ("alleles_%s.fa" % a.sample)
    out = a.work / "paf"
    out.mkdir(exist_ok=True)
    used = []
    for r, v in verified_runs(a):
        preset = PRESET[r["read_type"]]
        with open(str(out / ("%s_%s.paf" % (a.sample, r["run_accession"]))), "w") as h:
            subprocess.run(["minimap2", "-x", preset] + MM + ["-t", str(a.threads), str(fa), v["local_file"]], stdout=h,
                           stderr=subprocess.DEVNULL, check=True)
        used.append(dict(run=r["run_accession"], read_type=r["read_type"], preset=preset, sha256=v["sha256"]))
    (out / ("%s.json" % a.sample)).write_text(json.dumps(dict(sample=a.sample, runs=used, minimap2_args=MM), indent=1) + "\n")


def record_ids(header, grouping):
    """Molecule identifiers of one FASTQ header (without '@'): one key, or two for an ONT duplex record."""
    f = header.split()
    name, comment = f[0], (f[1] if len(f) > 1 else "")
    if grouping == "zmw_from_read_name":
        return ["/".join(comment.split("/")[:2])]
    if ONT_ID.match(comment):
        return comment.split("/")[0].split(";")  # ONT read id(s); a duplex "id1;id2" shares ids with its simplex parents
    return [name]  # SRA spot (subreads of one spot share the name) or no molecule information


def molecules(records):
    """records: [(read name, [ids])] -> {read name: molecule key}; records sharing any id form one molecule."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for _, ids in records:
        for i in ids[1:]:
            parent[find(i)] = find(ids[0])
    return {n: find(ids[0]) for n, ids in records}


def cmd_genome(a):
    fa = dict(x.split("=", 1) for x in a.fasta)[a.sample]
    out = a.work / "paf"
    used = []
    for r, v in verified_runs(a):
        run = r["run_accession"]
        allele_paf = out / ("%s_%s.paf" % (a.sample, run))
        if not allele_paf.exists():
            continue
        names = {line.split("\t", 1)[0] for line in open(str(allele_paf))}
        grouping = v["summary"]["molecule_grouping"]
        sub = a.work / ("subset_%s_%s.fq" % (a.sample, run))
        recs, keep = [], False
        with gzip.open(v["local_file"], "rt") as h, open(str(sub), "w") as o:
            for i, line in enumerate(h):
                if i % 4 == 0:
                    name = line[1:].split(None, 1)[0]
                    keep = name in names
                    if keep:
                        recs.append((name, record_ids(line[1:].strip(), grouping)))
                if keep:
                    o.write(line)
        mol = molecules(recs)
        write(out / ("molecules_%s_%s.tsv" % (a.sample, run)), [dict(read=n, molecule=m) for n, m in sorted(mol.items())])
        preset = PRESET[r["read_type"]]
        with open(str(out / ("genome_%s_%s.paf" % (a.sample, run))), "w") as h:
            subprocess.run(["minimap2", "-x", preset] + MM_GENOME + ["-t", str(a.threads), fa, str(sub)], stdout=h,
                           stderr=subprocess.DEVNULL, check=True)
        sub.unlink()
        used.append(dict(run=run, molecule_grouping=grouping, read_records=len(recs), read_names=len(mol),
                         molecules=len(set(mol.values())),
                         count_unit="read_records" if grouping == "molecule_grouping_unavailable" else "molecules"))
    (out / ("genome_%s.json" % a.sample)).write_text(json.dumps(dict(sample=a.sample, assembly=fa, runs=used,
                                                                     minimap2_args=MM_GENOME), indent=1) + "\n")


def parse_paf(path):
    for line in open(str(path)):
        f = line.rstrip("\n").split("\t")
        tags = dict(x.split(":", 2)[::2] for x in f[12:])
        yield dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], ts=int(f[7]), te=int(f[8]),
                   score=int(tags.get("AS", 0)), tp=tags.get("tp", "P"), cg=tags.get("cg", ""))


def cigar_ok(b, lo, hi, limit):
    """CIGAR check over the tested interval [lo, hi): merged indel events and the net indel excess must stay < limit.
    The read's own net indel rate outside the interval is subtracted (CLR/ONT reads have a systematic I/D bias)."""
    r = b["ts"]
    events, cur = [], None  # cur = [inserted, deleted] of the open event inside the interval
    net_in = net_out = out_len = 0
    for n, o in CIGAR.findall(b["cg"]):
        n = int(n)
        if o in "M=X":
            if cur and n > MERGE_MATCH:
                events.append(cur)
                cur = None
            out_len += n - overlap(r, r + n, lo, hi)
            r += n
            continue
        ins, dl = (n, 0) if o == "I" else (0, n) if o == "D" else (0, 0)
        if ins and lo <= r <= hi or dl and overlap(r, r + dl, lo, hi):
            net_in += ins - dl
            cur = [cur[0] + ins, cur[1] + dl] if cur else [ins, dl]
        else:
            net_out += ins - dl
            out_len += dl
            if cur:
                events.append(cur)
                cur = None
        r += dl
    if cur:
        events.append(cur)
    rate = net_out / out_len if out_len >= 2 * ANCHOR else 0.0
    largest = max([max(e) for e in events] or [0])
    return largest < limit and abs(net_in - rate * (hi - lo)) < limit


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def q_at(b, rpos):
    """Forward read coordinate of the base aligned at reference position rpos (deletion: next read base)."""
    r, o = b["ts"], 0
    for n, op in CIGAR.findall(b["cg"]):
        n = int(n)
        if op in "M=X":
            if r <= rpos < r + n:
                o += rpos - r
                break
            r, o = r + n, o + n
        elif op == "D":
            if r <= rpos < r + n:
                break
            r += n
        elif op == "I":
            o += n
    return b["qs"] + o if b["strand"] == "+" else b["qe"] - 1 - o


def r_at(g, q):
    """Reference position of forward read coordinate q in alignment g (insertion: the following reference base)."""
    o = q - g["qs"] if g["strand"] == "+" else g["qe"] - 1 - q
    r, qo = g["ts"], 0
    for n, op in CIGAR.findall(g["cg"]):
        n = int(n)
        if op in "M=X":
            if qo <= o < qo + n:
                return r + o - qo
            r, qo = r + n, qo + n
        elif op == "I":
            if qo <= o < qo + n:
                return r
            qo += n
        elif op == "D":
            r += n
    return None


def locus_check(alns, contig, strand, anchors):
    """ok / off_locus / locus_ambiguous for one read: anchors = [(read position, expected assembly position)]."""
    for q, exp in anchors:
        place = [(g, r_at(g, q)) for g in alns if g["qs"] <= q < g["qe"]]
        good = [g for g, r in place if g["t"] == contig and g["strand"] == strand and r is not None and abs(r - exp) <= LOCUS_TOL]
        if not good:
            return "off_locus"
        h = max(g["score"] for g in good)
        if any(all(g is not x for x in good) and g["score"] >= h - MIN_DELTA_AS for g, _ in place):
            return "locus_ambiguous"
    return "ok"


def test_geometry(t):
    """Per role: tested interval (lo, hi, middle length) per allele, CIGAR limit, read length needed for a
    length-matched count, and per allele the contig and expected assembly positions of the two anchors (lo, hi - 1):
    U side s - j1 + lo, D side e + (hi - 1 - j2), with (s, e) the site in the sample's own assembly."""
    geo = {}
    for role in "CE":
        contig, se = t["%s_site" % role].rsplit(":", 1)
        s, e = map(int, se.split("-"))
        iv, expect = {}, {}
        for al in "PA":
            j1, j2, sl, sr = map(int, t["%s%s_j" % (role, al)].split(","))
            lo, hi = max(0, j1 - ANCHOR - sl), min(int(t["%s%s_len" % (role, al)]), j2 + ANCHOR + sr)
            iv[al] = (lo, hi, j2 - j1)
            expect[al] = (s - j1 + lo, e + (hi - 1 - j2))
        geo[role] = dict(iv=iv, limit=max(MIN_INDEL, abs(iv["P"][2] - iv["A"][2]) / 2.0),
                         need=max(iv["P"][1] - iv["P"][0], iv["A"][1] - iv["A"][0]), contig=contig, expect=expect)
    return geo


def read_calls(best, galns, geo, qlen):
    """One read record -> {(test_id, role): call}; call is P/A (+ '_matched'), ambiguous, cigar_fail, off_locus or
    locus_ambiguous. best[(tid, role, allele)] = (2 spans and CIGAR ok / 1 spans, CIGAR fails / 0 partial, score,
    (strand, read position of anchor lo, of anchor hi - 1) when 2)."""
    out = {}
    for tid, role in {(k[0], k[1]) for k in best}:
        sp, sa = best.get((tid, role, "P"), (0, 0, None)), best.get((tid, role, "A"), (0, 0, None))
        if sp[0] < 2 and sa[0] < 2:
            if max(sp[0], sa[0]) == 1:
                out[(tid, role)] = "cigar_fail"
            continue
        if sp[0] == 2 and sa[0] == 2 and abs(sp[1] - sa[1]) < MIN_DELTA_AS:
            c = "ambiguous"
        elif sp[0] == 2 and (sa[0] < 2 or sp[1] > sa[1]):
            c = "P"
        else:
            c = "A"
        if galns is not None:
            al = "P" if sp[0] == 2 and (sa[0] < 2 or sp[1] >= sa[1]) else "A"  # anchors of the better-scoring allele
            strand, ql, qr = (sp if al == "P" else sa)[2]
            el, er = geo[tid][role]["expect"][al]
            lc = locus_check(galns, geo[tid][role]["contig"], strand, [(ql, el), (qr, er)])
            c = c if lc == "ok" else lc
        if c in ("P", "A") and qlen >= geo[tid][role]["need"]:
            c += "_matched"
        out[(tid, role)] = c
    return out


def classify(tests, runs):
    """runs: [(allele PAF, genome PAF or None, {read: molecule} or None)]. Returns (counts[(test_id, role)], number
    of molecules supporting 0 / 1 / >= 2 loci)."""
    geo = {t["test_id"]: test_geometry(t) for t in tests}
    test_locus = {t["test_id"]: t["locus_id"] for t in tests}
    counts = defaultdict(Counter)
    multi = Counter()
    for allele_paf, genome_paf, mol in runs:
        best, qlen = defaultdict(dict), {}
        for b in parse_paf(allele_paf):
            tid, role, al = b["t"].split("|")
            lo, hi, _ = geo[tid][role]["iv"][al]
            st = (2 if cigar_ok(b, lo, hi, geo[tid][role]["limit"]) else 1) if b["ts"] <= lo and b["te"] >= hi else 0
            qlen[b["q"]] = b["qlen"]
            cur = best[b["q"]].get((tid, role, al))
            if cur is None or (st, b["score"]) > cur[:2]:
                best[b["q"]][(tid, role, al)] = (st, b["score"], (b["strand"], q_at(b, lo), q_at(b, hi - 1)) if st == 2 else None)
        galns = defaultdict(list)
        if genome_paf:
            for g in parse_paf(genome_paf):
                galns[g["q"]].append(g)
        per_mol = defaultdict(lambda: defaultdict(set))  # molecule -> (tid, role) -> read-level calls
        for name, hits in best.items():
            m = mol.get(name, name) if mol else name
            for k, c in read_calls(hits, galns[name] if genome_paf else None, geo, qlen[name]).items():
                per_mol[m][k].add(c)
        for per in per_mol.values():
            loci = set()
            for k, cs in per.items():
                base = {c.split("_matched")[0] for c in cs}
                if len(base & {"P", "A", "ambiguous"}) > 1:
                    counts[k]["ambiguous"] += 1  # records of one molecule disagree
                    continue
                c = next(x for x in ("P", "A", "ambiguous", "locus_ambiguous", "off_locus", "cigar_fail") if x in base)
                counts[k][c] += 1
                if c in ("P", "A"):
                    loci.add(test_locus[k[0]])
                    if c + "_matched" in cs:
                        counts[k][c + "_matched"] += 1
            multi[min(len(loci), 2)] += 1
    return counts, multi


def verdict(expected, other):
    n = expected + other
    if n < MIN_SPAN:
        return "insufficient"
    if other <= MAX_OTHER_FRACTION * n:
        return "consistent"
    if expected <= MAX_OTHER_FRACTION * n:
        return "contradicted"
    return "mixed"


def test_verdict(cv, ev):
    if "contradicted" in (cv, ev):
        return "contradicted"
    if "mixed" in (cv, ev):
        return "mixed"
    if cv == ev == "consistent":
        return "both_samples_consistent"
    if cv == "consistent":
        return "carrier_only_consistent"
    if ev == "consistent":
        return "partner_only_consistent"
    return "not_tested"


def evidence_level(vs):
    if vs["contradicted"] or vs["mixed"]:
        return "X"
    if vs["both_samples_consistent"]:
        return "E2" if vs["both_samples_consistent"] == sum(vs.values()) else "E2s"
    if vs["carrier_only_consistent"] and vs["partner_only_consistent"]:
        return "E1cp"
    return "E1c" if vs["carrier_only_consistent"] else "E1p" if vs["partner_only_consistent"] else "E0"


def cmd_summarize(a):
    tests = read(a.work / "tests.tsv")
    paf_dir = a.work / "paf"
    have, unit, runs = set(), {}, []
    for p in sorted(paf_dir.glob("genome_*.json")):
        d = json.load(open(str(p)))
        if not d["runs"]:
            continue
        have.add(d["sample"])
        unit[d["sample"]] = "read_records" if any(r["count_unit"] == "read_records" for r in d["runs"]) else "molecules"
        for r in d["runs"]:
            mol = {x["read"]: x["molecule"] for x in read(paf_dir / ("molecules_%s_%s.tsv" % (d["sample"], r["run"])))}
            runs.append((paf_dir / ("%s_%s.paf" % (d["sample"], r["run"])), paf_dir / ("genome_%s_%s.paf" % (d["sample"], r["run"])), mol))
    counts, multi = classify(tests, runs)
    rows = []
    for t in tests:
        r = dict(t)
        for role, sample, exp in (("C", t["carrier"], "P"), ("E", t["empty"], "A")):
            c = counts[(t["test_id"], role)]
            oth = "A" if exp == "P" else "P"
            r.update({"%s_reads" % role: "yes" if sample in have else "no", "%s_count_unit" % role: unit.get(sample, ""),
                      "%s_present" % role: c["P"], "%s_absent" % role: c["A"],
                      "%s_present_matched" % role: c["P_matched"], "%s_absent_matched" % role: c["A_matched"],
                      "%s_ambiguous" % role: c["ambiguous"], "%s_locus_ambiguous" % role: c["locus_ambiguous"],
                      "%s_off_locus" % role: c["off_locus"], "%s_cigar_fail" % role: c["cigar_fail"],
                      "%s_verdict" % role: verdict(c[exp + "_matched"], c[oth + "_matched"]) if sample in have else "no_reads"})
        r["test_verdict"] = test_verdict(r["C_verdict"], r["E_verdict"])
        rows.append({k: v for k, v in r.items() if not k.endswith("_j")})
    write(a.out / "readcheck_tests.tsv", rows)
    by = defaultdict(list)
    for r in rows:
        by[r["locus_id"]].append(r)
    loci = []
    for lid, rs in by.items():
        vs = Counter(r["test_verdict"] for r in rs)
        loci.append(dict(locus_id=lid, carrier=rs[0]["carrier"], core_len=rs[0]["core_len"], locus_support=rs[0]["locus_support"],
                         boundary_consistency=rs[0]["boundary_consistency"],
                         difference_class=max((r["difference_class"] for r in rs), key=DIFF_ORDER.index),
                         tests=len(rs), test_verdicts=json.dumps(dict(vs)), evidence_level=evidence_level(vs),
                         carrier_count_unit=rs[0]["C_count_unit"],
                         partner_count_units=";".join(sorted({r["E_count_unit"] for r in rs if r["E_count_unit"]}))))
    write(a.out / "readcheck_loci.tsv", loci)
    short = [l for l in loci if int(l["core_len"]) < 2500]

    def tab(ls, key):
        return {s: dict(Counter(l["evidence_level"] for l in ls if l[key] == s)) for s in sorted({l[key] for l in ls})}

    excluded = Counter()
    for c in counts.values():
        excluded.update({k: v for k, v in c.items() if k in ("ambiguous", "locus_ambiguous", "off_locus", "cigar_fail")})
    summ = dict(caller=VERSION, flank=FLANK, anchor=ANCHOR, min_delta_as=MIN_DELTA_AS, min_span=MIN_SPAN,
                max_other_fraction=MAX_OTHER_FRACTION, min_indel=MIN_INDEL, locus_tol=LOCUS_TOL,
                tests=len(rows), test_verdicts=dict(Counter(r["test_verdict"] for r in rows)),
                side_verdicts=dict(Counter("%s/%s" % (r["C_verdict"], r["E_verdict"]) for r in rows)),
                loci=len(loci), evidence_levels=dict(Counter(l["evidence_level"] for l in loci)),
                short_loci=len(short), short_evidence_levels=dict(Counter(l["evidence_level"] for l in short)),
                short_by_difference_class=tab(short, "difference_class"), short_by_support=tab(short, "locus_support"),
                short_by_boundary=tab(short, "boundary_consistency"), short_by_carrier_count_unit=tab(short, "carrier_count_unit"),
                molecules_by_loci_supported={"0": multi[0], "1": multi[1], ">=2": multi[2]},
                excluded_molecule_test_pairs=dict(excluded), samples_with_reads=sorted(have), count_units=unit)
    (a.out / "readcheck_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["alleles", "fetch", "map", "genome", "summarize"])
    p.add_argument("--loci", type=Path)
    p.add_argument("--calls", type=Path)
    p.add_argument("--fasta", nargs="*")
    p.add_argument("--pool-dir", type=Path)
    p.add_argument("--reads", type=Path)
    p.add_argument("--work", type=Path)
    p.add_argument("--sample")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--jobs", type=int, default=1, help="fetch: concurrent run downloads")
    p.add_argument("--loci-ids", type=Path, help="alleles: TSV whose first column (after a header) lists the loci to test")
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    dict(alleles=cmd_alleles, fetch=cmd_fetch, map=cmd_map, genome=cmd_genome, summarize=cmd_summarize)[a.mode](a)


if __name__ == "__main__":
    main()
