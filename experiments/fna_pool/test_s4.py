"""Challenge cases for s4_annotate.py GFF3 handling (no external tools). Run: python test_s4.py (or pytest)."""

import random
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s4_annotate as s4  # noqa: E402
from s3_pair_calls import rc  # noqa: E402

R = random.Random(4)
SENSE = [c for c in s4.CODON if s4.CODON[c] != "*"]


def orf(n_codons):
    return "ATG" + "".join(R.choice(SENSE) for _ in range(n_codons - 2)) + "TAA"


def check(name, cond, detail):
    print(("PASS" if cond else "FAIL"), name, detail)
    return bool(cond)


def cases():
    L = 26000
    g = list("".join(R.choice("ACGT") for _ in range(L)))
    gff, faa = [], {}

    def put(pos, s):
        g[pos:pos + len(s)] = list(s)

    def line(uid, start, end, strand, phase=0, extra="", tt="11"):
        tt = ";transl_table=" + tt if tt else ""
        gff.append("c\tRefSeq\tCDS\t%d\t%d\t.\t%s\t%d\tID=%s;%s%s" % (start + 1, end, strand, phase, uid, extra, tt))

    a = orf(100)  # plus strand
    put(1000, a)
    line("cds-A", 1000, 1300, "+", 0, "protein_id=PA;product=plain")
    faa["PA"] = s4.translate(a)[:-1]
    b = orf(80)  # minus strand
    put(3000, rc(b))
    line("cds-B", 3000, 3240, "-", 0, "protein_id=PB;product=minus")
    faa["PB"] = s4.translate(b)[:-1]
    t = list(orf(100))  # -1 ribosomal frameshift: genome lacks one copy of base 149
    t[150:153] = list(t[149] + "CC")
    t = "".join(t)
    put(5000, t[:150] + t[151:])
    for s, e in ((5000, 5150), (5149, 5299)):
        line("cds-C", s, e, "+", 0, "protein_id=PC;exception=ribosomal slippage;product=frameshift")
    faa["PC"] = s4.translate(t)[:-1]
    d = "G" + orf(60)  # 5' partial, phase 1
    put(7000, d)
    line("cds-D", 7000, 7000 + len(d), "+", 1, "protein_id=PD;partial=true;start_range=.,7001;product=partial")
    faa["PD"] = s4.translate(d[1:])[:-1]
    o = orf(100)  # crosses the origin of a circular contig
    put(L - 150, o[:150])
    put(0, o[150:])
    line("cds-E", L - 150, L + 150, "+", 0, "protein_id=PE;product=origin")
    faa["PE"] = s4.translate(o)[:-1]
    put(9000, orf(50))
    line("cds-F", 9000, 9150, "+", 0, "pseudo=true;product=pseudo")
    sec = list(orf(60))
    sec[90:93] = list("TGA")
    sec = "".join(sec)
    put(11000, sec)
    line("cds-G", 11000, 11180, "+", 0, "protein_id=PG;transl_except=(pos:11091..11093,aa:Sec);product=selenoprotein")
    faa["PG"] = s4.translate(sec)[:-1].replace("*", "U")
    sec2 = list(orf(60))  # review case: two internal TGA, only the first is a stated exception
    sec2[90:93] = list("TGA")
    sec2[120:123] = list("TGA")
    sec2 = "".join(sec2)
    put(17000, sec2)
    line("cds-G2", 17000, 17180, "+", 0, "protein_id=PG2;transl_except=(pos:17091..17093,aa:Sec);product=one exception")
    faa["PG2"] = s4.translate(sec2)[:-1].replace("*", "U")
    sec3 = list(orf(60))  # minus strand exception
    sec3[60:63] = list("TGA")
    sec3 = "".join(sec3)
    put(18000, rc(sec3))
    line("cds-G3", 18000, 18180, "-", 0, "protein_id=PG3;transl_except=(pos:complement(18118..18120),aa:Sec);product=minus sec")
    faa["PG3"] = s4.translate(sec3)[:-1].replace("*", "U")
    h = orf(70)
    put(13000, h)
    line("cds-H", 13000, 13210, "+", 0, "protein_id=PH;Parent=gene-x,gene-y;product=multi parent")
    faa["PH"] = s4.translate(h)[:-1]
    m = orf(70)
    put(15000, m)
    line("cds-M", 15000, 15210, "+", 0, "protein_id=PM;product=faa differs")
    faa["PM"] = "M" + "W" * 68
    t4 = list(orf(80))  # Mycoplasma-style table 4: TTA start, internal TGA = Trp
    t4[0:3] = list("TTA")
    t4[60:63] = list("TGA")
    t4[150:153] = list("TGA")
    t4 = "".join(t4)
    p4 = "M" + "".join(s4.CODES[4]["codon"][t4[i:i + 3]] for i in range(3, len(t4) - 3, 3))
    put(20000, t4)
    line("cds-T4", 20000, 20240, "+", 0, "protein_id=P4;product=table 4", tt="4")
    faa["P4"] = p4
    put(20500, t4)
    line("cds-T4as11", 20500, 20740, "+", 0, "protein_id=P4b;product=table 4 gene labelled 11", tt="11")
    faa["P4b"] = p4
    put(21000, orf(60))
    line("cds-NOTT", 21000, 21180, "+", 0, "protein_id=PN;product=no table", tt="")
    faa["PN"] = s4.translate("".join(g[21000:21180]))[:-1]
    put(21500, orf(60))
    line("cds-T25", 21500, 21680, "+", 0, "protein_id=P25;product=unsupported table", tt="25")
    faa["P25"] = s4.translate("".join(g[21500:21680]))[:-1]
    put(22000, orf(60))
    line("cds-CONF", 22000, 22090, "+", 0, "protein_id=PX;product=conflict", tt="11")
    line("cds-CONF", 22090, 22180, "+", 0, "protein_id=PX;product=conflict", tt="4")
    faa["PX"] = s4.translate("".join(g[22000:22180]))[:-1]
    seq = "".join(g)
    tmp = Path(tempfile.mkdtemp(prefix="s4test_"))
    (tmp / "g.gff").write_text("##gff-version 3\n" + "\n".join(gff) + "\n")
    units = {u["id"]: u for u in s4.parse_gff(tmp / "g.gff", {"c": L})["c"]}
    res = {k: s4.coding_check(u, seq, faa) for k, u in units.items()}
    ok = []
    ok.append(check("plus and minus strand translate to the FAA protein", res["cds-A"] == ("match", []) and res["cds-B"] == ("match", []),
                    (res["cds-A"], res["cds-B"])))
    ok.append(check("multi-segment frameshift CDS is one unit and translates", len(units["cds-C"]["segments"]) == 2
                    and res["cds-C"] == ("match", []) and units["cds-C"]["exception"] == "ribosomal slippage", res["cds-C"]))
    ok.append(check("5' partial with phase 1: partial flag, translation matches", units["cds-D"]["partial"] and res["cds-D"][0] == "match",
                    res["cds-D"]))
    ok.append(check("origin-crossing unit: flagged, two pieces, translates", units["cds-E"]["crosses_origin"]
                    and units["cds-E"]["intervals"] == [(L - 150, L), (0, 150)] and res["cds-E"] == ("match", []),
                    (units["cds-E"]["intervals"], res["cds-E"])))
    ok.append(check("pseudo unit is not compared", res["cds-F"][0] == "pseudo", res["cds-F"]))
    ok.append(check("transl_except selenocysteine", res["cds-G"] == ("match_with_transl_except", []), res["cds-G"]))
    ok.append(check("unstated second internal stop is not excused by transl_except", res["cds-G2"][0] == "mismatch"
                    and "internal_stop" in res["cds-G2"][1], res["cds-G2"]))
    ok.append(check("minus-strand transl_except placed on its codon", res["cds-G3"] == ("match_with_transl_except", []), res["cds-G3"]))
    ok.append(check("multi-Parent kept", units["cds-H"]["parent"] == "gene-x,gene-y" and res["cds-H"][0] == "match", units["cds-H"]["parent"]))
    ok.append(check("FAA differs -> mismatch", res["cds-M"][0] == "mismatch" and "translation_mismatch" in res["cds-M"][1], res["cds-M"]))
    ok.append(check("table 4 (TTA start, TGA=Trp) translates under its own transl_table", res["cds-T4"] == ("match", []), res["cds-T4"]))
    ok.append(check("table 4 gene labelled 11 is not accepted", res["cds-T4as11"][0] == "mismatch"
                    and "internal_stop" in res["cds-T4as11"][1] and "no_start" in res["cds-T4as11"][1], res["cds-T4as11"]))
    ok.append(check("missing, unsupported or conflicting transl_table -> review, not translated",
                    res["cds-NOTT"] == ("transl_table_review", ["transl_table_missing"])
                    and res["cds-T25"] == ("transl_table_review", ["unsupported_transl_table=25"])
                    and res["cds-CONF"] == ("transl_table_review", ["transl_table_conflict=11/4"]),
                    (res["cds-NOTT"], res["cds-T25"], res["cds-CONF"])))
    us = list(units.values())
    rel = {k: sorted(r for r, _ in s4.boundary_relations(us, s, e)) for k, (s, e) in
           dict(inside_gene=(1100, 1200), most_of_gene=(1050, 1400), edge_of_gene=(1250, 1600), clear=(1400, 1600),
                origin_piece=(100, 400)).items()}
    ok.append(check("boundary relations", rel == dict(inside_gene=["host_gene_interrupted"], most_of_gene=["gene_mostly_inside"],
                                                      edge_of_gene=["gene_mostly_outside"], clear=[], origin_piece=["gene_mostly_outside"]), rel))
    return ok


def test_all():
    assert all(cases())


if __name__ == "__main__":
    r = cases()
    print("%d/%d passed" % (sum(r), len(r)))
    sys.exit(0 if all(r) else 1)
