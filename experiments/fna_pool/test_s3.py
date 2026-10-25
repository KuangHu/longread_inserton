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
    return ok


def test_all():
    assert all(bool(x) for x in cases())


if __name__ == "__main__":
    res = cases()
    print("%d/%d passed" % (sum(bool(x) for x in res), len(res)))
    sys.exit(0 if all(bool(x) for x in res) else 1)
