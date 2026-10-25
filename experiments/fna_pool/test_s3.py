"""Synthetic challenge cases for s3_pair_calls.py (needs minimap2 on PATH). Run: python test_s3.py  (or pytest)."""

import random
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s3_pair_calls as s3  # noqa: E402

R = random.Random(20261024)
BASE = "".join(R.choice("ACGT") for _ in range(300000))


def rand(n):
    return "".join(R.choice("ACGT") for _ in range(n))


def snps(seq, rate, seed):
    r = random.Random(seed)
    s = list(seq)
    for i in range(len(s)):
        if r.random() < rate:
            s[i] = r.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def ins(seq, pos, x):
    return seq[:pos] + x + seq[pos:]


def run(tmp, genomes, pairs):
    """genomes: name -> sequence (one contig 'c'); pairs: list of (a, b). Returns (calls, loci)."""
    paths = {}
    for n, s in genomes.items():
        paths[n] = tmp / (n + ".fa")
        paths[n].write_text(">c\n" + s + "\n")
    calls, one = [], []
    for i, (x, y) in enumerate(pairs):
        for ref, qry in ((x, y), (y, x)):
            paf = tmp / ("%s_on_%s.paf" % (qry, ref))
            s3.run_minimap2(paths[ref], paths[qry], paf, 2)
            c, o = s3.call_pair("T%d" % i, ref, {"c": genomes[ref]}, qry, {"c": genomes[qry]}, paf,
                                {ref: paths[ref], qry: paths[qry]}, tmp, 2)
            calls, one = calls + c, one + o
    calls = [{k: str(v) for k, v in c.items()} for c in calls]
    run.one_flank = one
    return calls, s3.merge_calls(calls)


def check(name, cond, detail):
    print(("PASS" if cond else "FAIL"), name, detail)
    return bool(cond)


def cases():
    tmp = Path(tempfile.mkdtemp(prefix="s3test_"))
    B = snps(BASE, 0.005, 1)
    X = rand(1500)
    ok = []

    calls, loci = run(tmp, dict(A=ins(BASE, 100000, X), B=B), [("A", "B")])
    L = [l for l in loci if l["carrier"] == "A"]
    ok.append(check("simple insertion: one anchored locus, both mappings", len(loci) == 1 and L and L[0]["status"] == "two_flank_anchored"
                    and abs(L[0]["core_start"] - 100000) <= 2 and L[0]["core_len"] == 1500 and L[0]["mappings"] == 2,
                    [(l["carrier"], l["core_start"], l["core_len"], l["status"], l["mappings"]) for l in loci]))

    calls, loci = run(tmp, dict(A=ins(BASE, 100000, X), B=s3.rc(B)), [("A", "B")])
    ok.append(check("reverse-complemented partner", len(loci) == 1 and loci[0]["status"] == "two_flank_anchored"
                    and abs(loci[0]["core_start"] - 100000) <= 2 and {c["aln_strand"] for c in calls} == {"-"},
                    [(l["carrier"], l["core_start"], l["core_len"], l["status"]) for l in loci]))

    # Read-check alleles on the same minus-strand pair, from both mapping directions: the partner's alleles must be
    # the reverse complement of the carrier's (flanks are identical apart from 0.5% SNPs, so compare against BASE).
    import s3_readcheck as rcheck
    Aseq, Bseq = ins(BASE, 100000, X), s3.rc(BASE)
    calls, _ = run(tmp, dict(A=Aseq, B=Bseq), [("A", "B")])
    res = []
    for c in calls:
        al = rcheck.build_alleles(c, Aseq, Bseq)
        cp, ca, ep, ea = al[("C", "P")], al[("C", "A")], al[("E", "P")], al[("E", "A")]
        res.append(cp[0] == s3.rc(ep[0]) and ca[0] == s3.rc(ea[0]) and ca[0] in BASE
                   and cp[2] - cp[1] == 1500 and ea[2] == ea[1] and ep[2] - ep[1] == 1500)
    ok.append(check("read-check alleles, minus-strand partner", len(res) == 2 and all(res), res))

    X2 = rand(700)
    A = ins(ins(BASE, 100000, X), 100000 + 1500 + 300, X2)
    calls, loci = run(tmp, dict(A=A, B=B), [("A", "B")])
    ok.append(check("adjacent insertions 300 bp apart stay separate", sorted(l["core_len"] for l in loci) == [700, 1500]
                    and all(l["status"] == "two_flank_anchored" for l in loci),
                    [(l["core_start"], l["core_len"], l["status"]) for l in loci]))

    T = BASE[100000 - 800:100000]
    calls, loci = run(tmp, dict(A=ins(BASE, 100000, T), B=B), [("A", "B")])
    ok.append(check("tandem duplication flagged", len(loci) == 1 and loci[0]["tandem_like"] and loci[0]["core_len"] == 800,
                    [(l["core_start"], l["core_len"], l["status"], l["tandem_like"]) for l in loci]))

    Bd = ins(B, 200000, BASE[99500:100000])
    calls, loci = run(tmp, dict(A=ins(BASE, 100000, X), B=Bd), [("A", "B")])
    flagged = [c for c in calls if c["carrier"] == "A" and ("flank_not_unique" in c["reasons"] or "flank_competing" in c["reasons"])]
    ok.append(check("duplicated left flank -> not unique, not anchored", flagged and all(l["status"] != "two_flank_anchored"
                    for l in loci if l["carrier"] == "A" and l["core_len"] == 1500),
                    [(l["carrier"], l["core_start"], l["core_len"], l["status"], l["reasons"]) for l in loci]))

    calls, loci = run(tmp, dict(A=ins(BASE, 500, X), B=B), [("A", "B")])
    edge = [o for o in run.one_flank if o["near_contig_end"]]
    ok.append(check("near contig end -> no anchored locus, one-flank record kept", edge and not any(l["status"] == "two_flank_anchored" for l in loci),
                    ([(l["core_start"], l["core_len"], l["status"]) for l in loci],
                     [(o["mapping"], o["side"], o["ref_unaligned_beyond"], o["qry_unaligned_beyond"]) for o in run.one_flank])))

    Y = rand(800)
    calls, loci = run(tmp, dict(A=ins(BASE, 100000, X[:1000]), B=ins(B, 100000, Y)), [("A", "B")])
    ok.append(check("replacement X vs Y", loci and all(l["status"] == "replacement" for l in loci),
                    [(l["carrier"], l["core_start"], l["core_len"], l["status"]) for l in loci]))

    calls, loci = run(tmp, dict(A=ins(BASE, 100000, X), B=ins(B, 200000, X)), [("A", "B")])
    by = {l["carrier"]: l for l in loci}
    ok.append(check("same core at different sites -> two loci, core copy in the other genome", set(by) == {"A", "B"}
                    and all(int(l["core_copies_empty"]) >= 1 and l["status"] == "two_flank_anchored" for l in loci),
                    [(l["carrier"], l["core_start"], l["core_len"], l["status"], l["core_copies_carrier"], l["core_copies_empty"]) for l in loci]))

    C = snps(BASE, 0.005, 2)
    calls, loci = run(tmp, dict(A=ins(BASE, 100000, X), B=B, C=C), [("A", "B"), ("A", "C")])
    ok.append(check("one locus from two pairs merges", len(loci) == 1 and loci[0]["n_pairs"] == 2,
                    [(l["carrier"], l["core_start"], l["n_pairs"], l["n_calls"]) for l in loci]))

    calls, loci = run(tmp, dict(A=ins(ins(BASE, 150000, rand(30)), 100000, rand(60)), B=B), [("A", "B")])
    ok.append(check("60 bp indel called, 30 bp not", sorted(l["core_len"] for l in loci) == [60],
                    [(l["core_start"], l["core_len"], l["status"], l["n_calls"]) for l in loci]))

    big = rand(20000)
    calls, loci = run(tmp, dict(A=ins(BASE, 100000, big), B=B), [("A", "B")])
    ok.append(check("20 kb insertion (split alignment)", len(loci) == 1 and abs(loci[0]["core_start"] - 100000) <= 2
                    and abs(loci[0]["core_len"] - 20000) <= 4 and loci[0]["status"] == "two_flank_anchored",
                    [(l["core_start"], l["core_len"], l["status"], Counter(c["source"] for c in calls)) for l in loci]))

    huge = rand(60000)
    calls, loci = run(tmp, dict(A=ins(BASE, 100000, huge), B=B), [("A", "B")])
    src = Counter(c["source"] for c in calls)
    ok.append(check("60 kb insertion", len(loci) == 1 and abs(loci[0]["core_start"] - 100000) <= 2
                    and abs(loci[0]["core_len"] - 60000) <= 4 and loci[0]["status"] == "two_flank_anchored",
                    [(l["core_start"], l["core_len"], l["status"], dict(src)) for l in loci]))

    # minimap2 bridged both large insertions inside one CIGAR; exercise the split-block path with a two-block PAF.
    A = ins(BASE, 100000, huge)
    paf = tmp / "blocks.paf"
    tags = "tp:A:P\tAS:i:1\tcg:Z:%d="
    paf.write_text("c\t300000\t0\t100000\t+\tc\t360000\t0\t100000\t100000\t100000\t60\t" + tags % 100000 + "\n"
                   "c\t300000\t100000\t300000\t+\tc\t360000\t160000\t360000\t200000\t200000\t60\t" + tags % 200000 + "\n")
    c, _ = s3.call_pair("TB", "A", {"c": A}, "B", {"c": BASE}, paf)
    ok.append(check("split-block path (two-block PAF)", len(c) == 1 and c[0]["source"] == "split_blocks" and c[0]["carrier"] == "A"
                    and (c[0]["core_start"], c[0]["core_end"]) == (100000, 160000) and c[0]["status"] == "two_flank_anchored",
                    [(x["source"], x["carrier"], x["core_start"], x["core_end"], x["status"], x["reasons"]) for x in c]))
    ok += readcheck_cases(tmp)
    return ok


def readcheck_cases(tmp):
    """Read classification on hand-written PAF records (no minimap2), incl. the 0ed1ca4 review counterexamples."""
    import s3_readcheck as rcheck
    ok = []

    def test(tid, locus, cs, core):
        return {"test_id": tid, "locus_id": locus, "CP_len": 4000 + core, "CP_j": "2000,%d,0,0" % (2000 + core),
                "CA_len": 4000, "CA_j": "2000,2000,0,0", "EP_len": 4000 + core, "EP_j": "2000,%d,0,0" % (2000 + core),
                "EA_len": 4000, "EA_j": "2000,2000,0,0", "C_site": "c:%d-%d" % (cs, cs + core), "E_site": "e:%d-%d" % (cs, cs)}

    def paf(name, lines):
        p = tmp / name
        p.write_text("".join("%s\t%d\t%d\t%d\t+\t%s\t%d\t%d\t%d\t1\t1\t60\ttp:A:%s\tAS:i:%d\tcg:Z:%s\n" % l for l in lines))
        return p

    tests = [test("T1", "L1", 10000, 54), test("T2", "L2", 50000, 54), test("T3", "L3", 20000, 54),
             test("T4", "L4", 20400, 54), test("T5", "L5", 30000, 5000)]
    al = paf("al.paf", [
        # r1: alignment to the present allele deletes the whole 54 bp core; no absent-allele record (review case 1)
        ("r1", 4000, 0, 4000, "T1|C|P", 4054, 0, 4054, "P", 7800, "2000M54D2000M"),
        # r2: strong on T1, weak on T2 (review case 2)
        ("r2", 4054, 0, 4054, "T1|C|P", 4054, 0, 4054, "P", 8000, "4054M"),
        ("r2", 4054, 0, 4054, "T2|C|P", 4054, 0, 4054, "S", 3000, "4054M"),
        # r3: one molecule across two neighbouring loci, both truly present
        ("r3", 7000, 0, 4054, "T3|C|P", 4054, 0, 4054, "P", 8000, "4054M"),
        ("r3", 7000, 2900, 6954, "T4|C|P", 4054, 0, 4054, "P", 8000, "4054M"),
        # r4: absent at T5 (5 kb core) but too short to have spanned the present allele
        ("r4", 3000, 0, 3000, "T5|C|A", 4000, 500, 3500, "P", 6000, "3000M"),
        # r5a/r5b: two records of one molecule, both present at T1
        ("r5a", 4054, 0, 4054, "T1|C|P", 4054, 0, 4054, "P", 8000, "4054M"),
        ("r5b", 4054, 0, 4054, "T1|C|P", 4054, 0, 4054, "P", 8000, "4054M")])
    ge = paf("ge.paf", [
        ("r1", 4000, 0, 4000, "c", 10 ** 6, 8000, 12054, "P", 7800, "4000M"),
        ("r2", 4054, 0, 4054, "c", 10 ** 6, 8000, 12054, "P", 8000, "4054M"),
        ("r2", 4054, 0, 4054, "c", 10 ** 6, 48000, 52054, "S", 3000, "4054M"),
        ("r3", 7000, 0, 7000, "c", 10 ** 6, 18000, 24454, "P", 14000, "7000M"),
        ("r4", 3000, 0, 3000, "c", 10 ** 6, 28500, 31500, "P", 6000, "3000M"),
        ("r5a", 4054, 0, 4054, "c", 10 ** 6, 8000, 12054, "P", 8000, "4054M"),
        ("r5b", 4054, 0, 4054, "c", 10 ** 6, 8000, 12054, "P", 8000, "4054M")])
    mol = {"r1": "r1", "r2": "r2", "r3": "r3", "r4": "r4", "r5a": "m5", "r5b": "m5"}
    counts, multi = rcheck.classify(tests, [(al, ge, mol)])
    c = {t: counts[(t, "C")] for t in ("T1", "T2", "T3", "T4", "T5")}
    ok.append(check("readcheck: deleted core in CIGAR is not present", c["T1"]["cigar_fail"] == 1 and c["T1"]["P"] == 2,
                    dict(c["T1"])))
    ok.append(check("readcheck: weak second locus loses to genome-wide placement", c["T2"]["P"] == 0 and c["T2"]["locus_ambiguous"] == 1,
                    dict(c["T2"])))
    ok.append(check("readcheck: one read across two neighbouring loci supports both", c["T3"]["P"] == c["T4"]["P"] == 1 and multi[2] == 1,
                    (dict(c["T3"]), dict(c["T4"]), dict(multi))))
    ok.append(check("readcheck: short absent read is raw but not length-matched", c["T5"]["A"] == 1 and c["T5"]["A_matched"] == 0,
                    dict(c["T5"])))
    v01 = Counter()
    for b in rcheck.parse_paf(al):  # the v0.1 rule (span only) would have counted r1 present and r2 at T2
        tid, role, a_ = b["t"].split("|")
        lo, hi, _ = rcheck.test_geometry([t for t in tests if t["test_id"] == tid][0])[role]["iv"][a_]
        v01[tid] += b["ts"] <= lo and b["te"] >= hi and a_ == "P"
    ok.append(check("readcheck: the same records pass the v0.1 span-only rule", v01["T1"] == 4 and v01["T2"] == 1, dict(v01)))
    dup = rcheck.molecules([("a", rcheck.record_ids("a 0f1e2d3c-0000-4000-8000-000000000001 st:Z:x", "one_record_per_molecule")),
                            ("b", rcheck.record_ids("b 0f1e2d3c-0000-4000-8000-000000000002 st:Z:x", "one_record_per_molecule")),
                            ("d", rcheck.record_ids("d 0f1e2d3c-0000-4000-8000-000000000001;0f1e2d3c-0000-4000-8000-000000000002 st", "one_record_per_molecule")),
                            ("z", rcheck.record_ids("z m64_1/77/ccs/1", "zmw_from_read_name"))])
    ok.append(check("readcheck: ONT duplex record and both simplex parents are one molecule",
                    dup["a"] == dup["b"] == dup["d"] and dup["z"] == "m64_1/77", dup))
    core = rand(300)
    dc = [s3.difference_class(core, ""), s3.difference_class(core, core[:3] + core[-2:]), s3.difference_class(core, rand(8)),
          s3.difference_class(core, rand(60))]
    ok.append(check("difference_class: absorbed ends, small other, replacement",
                    [d[0] for d in dc] == ["presence_absence", "presence_absence", "complex_small_other", "replacement"], dc))
    return ok


def test_all():
    assert all(bool(x) for x in cases())


if __name__ == "__main__":
    res = cases()
    print("%d/%d passed" % (sum(bool(x) for x in res), len(res)))
    sys.exit(0 if all(bool(x) for x in res) else 1)
