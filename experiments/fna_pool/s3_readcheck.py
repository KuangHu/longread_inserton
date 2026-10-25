"""S3 long-read allele check: do a sample's own long reads span the allele its assembly shows at a call?

Unit = one test: a distinct raw call (carrier core coordinates + partner empty site) of a selected locus. Each
sample is tested against two local alleles built from its own assembly, so flank divergence between strains cannot
favour one allele:
  carrier (role C): present = U + core + D, absent = U + other + D   (other = partner's other-side segment, oriented)
  partner (role E): absent = U' + other + D', present = U' + core + D' (core oriented to the partner's strand)
with FLANK bp of U/D. A read supports an allele when its alignment to that allele covers both junctions plus ANCHOR bp
(plus the boundary shift window) on each side; a read spanning both alleles goes to the higher alignment score if
the difference is >= MIN_DELTA_AS, otherwise it is ambiguous. Partial reads are not counted.

Role verdicts (expected allele: C present, E absent), from spanning reads n = expected + other:
  insufficient (n < MIN_SPAN), consistent (other <= MAX_OTHER_FRACTION * n), contradicted (expected <= that
  fraction), mixed (both alleles). These are fixed in advance and not tuned on the pilot.
The reads are those deposited for the assembly's BioSample and were usually used to build it: this is an
assembly-consistency check (misassembly, collapsed repeats, within-sample allele mixtures), not independent
sequencing validation.

alleles   tests.tsv + alleles_<sample>.fa from loci/calls (loci with status two_flank_anchored; all their calls).
fetch     Download the usable runs of the tested samples (scripts/fetch_long_reads.py fetch logic, MD5 checked).
map       minimap2 (preset by read type) of one sample's runs against its allele FASTA.
summarize Per test and per locus counts and verdicts.
"""

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from s3_pair_calls import read_fasta, rc  # noqa: E402

VERSION = "s3_readcheck_v0.1"
FLANK = 2000
ANCHOR = 300
MIN_DELTA_AS = 20
MIN_SPAN = 3
MAX_OTHER_FRACTION = 0.1
PRESET = {"ont": "map-ont", "pacbio_hifi_ccs": "map-hifi", "pacbio_clr_subreads": "map-pb", "pacbio_unresolved": "map-pb"}
MM = ["-c", "--secondary=yes", "-N", "100", "-p", "0.3"]


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
    loci = {l["locus_id"]: l for l in read(a.loci)}
    fa = dict(x.split("=", 1) for x in a.fasta)
    seqs = {}
    tests, alleles = {}, defaultdict(list)
    for c in read(a.calls):
        l = loci[c["locus_id"]]
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
        al = build_alleles(c, seqs[c["carrier"]][c["carrier_contig"]], seqs[c["empty"]][c["empty_contig"]])
        t = dict(test_id=tid, locus_id=c["locus_id"], pair_id=c["pair_id"], carrier=c["carrier"], empty=c["empty"],
                 core_len=int(c["core_end"]) - int(c["core_start"]), other_side_gap=c["other_side_gap"],
                 aln_strand=c["aln_strand"], call_status_first=c["status"], locus_support=l.get("anchored_support", ""),
                 boundary_consistency=l["boundary_consistency"], calls=1, call_status=Counter([c["status"]]))
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
    tests = read(a.work / "tests.tsv")
    samples = {t["carrier"] for t in tests} | {t["empty"] for t in tests}
    for s, rs in sorted(sample_runs(a.pool_dir, samples).items()):
        for r in rs:
            if (a.reads / r["run_accession"] / "verify.json").exists():
                print(s, r["run_accession"], "already fetched")
                continue
            if ";" in r["fastq_ftp"]:
                print(s, r["run_accession"], "skipped: more than one FASTQ file", file=sys.stderr)
                continue
            run = dict(run=r["run_accession"], sample_id=s, url="https://" + r["fastq_ftp"], fastq_bytes=int(r["fastq_bytes"]),
                       fastq_md5=r["fastq_md5"], platform=r["instrument_platform"], archive_base_count=int(r["base_count"]),
                       archive_read_count=int(r["read_count"]))
            res = flr.fetch(run, a.reads)
            print(s, r["run_accession"], res["status"])


def cmd_map(a):
    runs = sample_runs(a.pool_dir, {a.sample})[a.sample]
    fa = a.work / ("alleles_%s.fa" % a.sample)
    out = a.work / "paf"
    out.mkdir(exist_ok=True)
    used = []
    for r in runs:
        v = a.reads / r["run_accession"] / "verify.json"
        if not v.exists():
            print(a.sample, r["run_accession"], "not fetched", file=sys.stderr)
            continue
        v = json.load(open(str(v)))
        if v["status"] != "verified":
            print(a.sample, r["run_accession"], "not verified:", v["checks"], file=sys.stderr)
            continue
        preset = PRESET[r["read_type"]]
        with open(str(out / ("%s_%s.paf" % (a.sample, r["run_accession"]))), "w") as h:
            subprocess.run(["minimap2", "-x", preset] + MM + ["-t", str(a.threads), str(fa), v["local_file"]], stdout=h,
                           stderr=subprocess.DEVNULL, check=True)
        used.append(dict(run=r["run_accession"], read_type=r["read_type"], preset=preset, sha256=v["sha256"]))
    (out / ("%s.json" % a.sample)).write_text(json.dumps(dict(sample=a.sample, runs=used, minimap2_args=MM), indent=1) + "\n")


def parse_paf(path):
    for line in open(str(path)):
        f = line.rstrip("\n").split("\t")
        tags = dict(x.split(":", 2)[::2] for x in f[12:])
        yield dict(q=f[0], t=f[5], ts=int(f[7]), te=int(f[8]), score=int(tags.get("AS", 0)), tp=tags.get("tp", "P"))


def classify(paf_files, tests):
    """counts[(test_id, role)] = Counter(P/A/ambiguous)."""
    jt = {}
    for t in tests:
        for role in "CE":
            for al in "PA":
                j1, j2, sl, sr = map(int, t["%s%s_j" % (role, al)].split(","))
                jt[(t["test_id"], role, al)] = (j1 - ANCHOR - sl, j2 + ANCHOR + sr)
    best = defaultdict(dict)  # (read, run) -> {(tid, role, allele): (spans, score)}
    for p in paf_files:
        for b in parse_paf(p):
            tid, role, al = b["t"].split("|")
            lo, hi = jt[(tid, role, al)]
            spans = b["ts"] <= max(0, lo) and b["te"] >= hi
            k = (tid, role, al)
            cur = best[(b["q"], p.name)].get(k)
            if cur is None or (spans, b["score"]) > cur:
                best[(b["q"], p.name)][k] = (spans, b["score"])
    counts = defaultdict(Counter)
    for hits in best.values():
        for tid, role in {(k[0], k[1]) for k in hits}:
            sp, sa = hits.get((tid, role, "P"), (False, 0)), hits.get((tid, role, "A"), (False, 0))
            if not sp[0] and not sa[0]:
                continue
            if sp[0] and sa[0] and abs(sp[1] - sa[1]) < MIN_DELTA_AS:
                counts[(tid, role)]["ambiguous"] += 1
            elif sp[0] and (not sa[0] or sp[1] > sa[1]):
                counts[(tid, role)]["P"] += 1
            else:
                counts[(tid, role)]["A"] += 1
    return counts


def verdict(expected, other):
    n = expected + other
    if n < MIN_SPAN:
        return "insufficient"
    if other <= MAX_OTHER_FRACTION * n:
        return "consistent"
    if expected <= MAX_OTHER_FRACTION * n:
        return "contradicted"
    return "mixed"


def cmd_summarize(a):
    tests = read(a.work / "tests.tsv")
    paf_dir = a.work / "paf"
    have = {p.stem for p in paf_dir.glob("*.json")}
    counts = classify(sorted(paf_dir.glob("*.paf")), tests)
    rows = []
    for t in tests:
        r = dict(t)
        for role, sample, exp in (("C", t["carrier"], "P"), ("E", t["empty"], "A")):
            c = counts[(t["test_id"], role)]
            oth = "A" if exp == "P" else "P"
            r.update({"%s_reads" % role: "yes" if sample in have else "no", "%s_present" % role: c["P"],
                      "%s_absent" % role: c["A"], "%s_ambiguous" % role: c["ambiguous"],
                      "%s_verdict" % role: verdict(c[exp], c[oth]) if sample in have else "no_reads"})
        v = (r["C_verdict"], r["E_verdict"])
        if v == ("consistent", "consistent"):
            r["test_verdict"] = "both_sides_consistent"
        elif "contradicted" in v:
            r["test_verdict"] = "contradicted"
        elif "mixed" in v:
            r["test_verdict"] = "mixed"
        elif "consistent" in v:
            r["test_verdict"] = "one_side_consistent"
        else:
            r["test_verdict"] = "not_tested"
        rows.append({k: v for k, v in r.items() if not k.endswith("_j")})
    write(a.out / "readcheck_tests.tsv", rows)
    by = defaultdict(list)
    for r in rows:
        by[r["locus_id"]].append(r)
    loci = []
    for lid, rs in by.items():
        vs = Counter(r["test_verdict"] for r in rs)
        loci.append(dict(locus_id=lid, carrier=rs[0]["carrier"], core_len=rs[0]["core_len"], locus_support=rs[0]["locus_support"],
                         boundary_consistency=rs[0]["boundary_consistency"], tests=len(rs), test_verdicts=json.dumps(dict(vs)),
                         locus_verdict=("contradicted" if vs["contradicted"] else "mixed" if vs["mixed"] else
                                        "both_sides_consistent" if vs["both_sides_consistent"] else
                                        "one_side_consistent" if vs["one_side_consistent"] else "not_tested")))
    write(a.out / "readcheck_loci.tsv", loci)
    short = [l for l in loci if int(l["core_len"]) < 2500]
    summ = dict(caller=VERSION, flank=FLANK, anchor=ANCHOR, min_delta_as=MIN_DELTA_AS, min_span=MIN_SPAN,
                max_other_fraction=MAX_OTHER_FRACTION, tests=len(rows), test_verdicts=dict(Counter(r["test_verdict"] for r in rows)),
                loci=len(loci), locus_verdicts=dict(Counter(l["locus_verdict"] for l in loci)),
                short_loci=len(short), short_locus_verdicts=dict(Counter(l["locus_verdict"] for l in short)),
                short_by_support={s: dict(Counter(l["locus_verdict"] for l in short if l["locus_support"] == s))
                                  for s in sorted({l["locus_support"] for l in short})},
                short_partner_dependent=dict(Counter(l["locus_verdict"] for l in short if l["boundary_consistency"] == "partner_dependent")),
                samples_with_reads=sorted(have))
    (a.out / "readcheck_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["alleles", "fetch", "map", "summarize"])
    p.add_argument("--loci", type=Path)
    p.add_argument("--calls", type=Path)
    p.add_argument("--fasta", nargs="*")
    p.add_argument("--pool-dir", type=Path)
    p.add_argument("--reads", type=Path)
    p.add_argument("--work", type=Path)
    p.add_argument("--sample")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    dict(alleles=cmd_alleles, fetch=cmd_fetch, map=cmd_map, summarize=cmd_summarize)[a.mode](a)


if __name__ == "__main__":
    main()
