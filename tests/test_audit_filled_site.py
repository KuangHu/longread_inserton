"""Challenge tests with known truth for scripts/audit_filled_site.py (protocol T0)."""

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_filled_site.py"
X = (5000, 6200)  # candidate on c1, length 1200
EDGE = (50, 1250)  # candidate whose left anchor would run off the contig


def paf(q, ts, te, cg, tp="P", strand="+", t="c1", AS=7000, mapq=60, qlen=None):
    ref = te - ts
    qlen = qlen or ref
    return "\t".join(map(str, [q, qlen, 0, qlen, strand, t, 20000, ts, te, ref, ref, mapq,
                               "tp:A:" + tp, "AS:i:%d" % AS, "cg:Z:" + cg]))


def full(q, off=0, **kw):
    """A clean full span; distinct reads get distinct `off` so they are not merged as duplicates."""
    return paf(q, 4000 - off, 7500, "%dM" % (3500 + off), **kw)


class AuditChallengeTests(unittest.TestCase):
    def run_case(self, lines, recheck=None, members=None, expect_rc=0):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "ref.fa").write_text(">c1\n" + "A" * 20000 + "\n>c2\n" + "A" * 20000 + "\n")
            rows = [dict(member=name, family="F", method="RS", accession="c1", candidate_start0=s,
                         candidate_end0=e, boundary_source="provisional_homology_median")
                    for name, (s, e) in (members or {"m1": X}).items()]
            (tmp / "members.json").write_text(json.dumps(rows))
            (tmp / "reads.paf").write_text("".join(l + "\n" for l in lines))
            cmd = [sys.executable, str(SCRIPT), "classify", "--members", str(tmp / "members.json"),
                   "--reference", str(tmp / "ref.fa"), "--paf", str(tmp / "reads.paf"),
                   "--sample", "test", "--out", str(tmp / "out")]
            if recheck is not None:
                (tmp / "recheck.paf").write_text("".join(l + "\n" for l in recheck))
                cmd += ["--recheck-paf", str(tmp / "recheck.paf")]
            result = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(result.returncode, expect_rc, result.stderr)
            if expect_rc:
                return None
            check = subprocess.run([sys.executable, str(SCRIPT.with_name("check_audited_evidence.py")),
                                    "--audit-dir", str(tmp / "out"), "--members", str(tmp / "members.json"),
                                    "--reference", str(tmp / "ref.fa"), "--out", str(tmp / "check.json")],
                                   capture_output=True, text=True)
            self.assertEqual(check.returncode, 0, check.stderr)
            with open(tmp / "out/member_audit.tsv") as h:
                audit = {r["member"]: r for r in csv.DictReader(h, delimiter="\t")}
            self.assertEqual(json.loads((tmp / "out/summary.json").read_text())["schema"], "audited_filled_site_v1")
            self.assertEqual(json.loads((tmp / "out/summary.json").read_text())["classifier_revision"], "nearby_competition_v2")
            return audit

    def cls(self, lines, **kw):
        return self.run_case(lines, **kw)["m1"]["audited_class"]

    def test_full_span_two_reads(self):
        self.assertEqual(self.cls([full("r1"), full("r2", off=7, AS=6990)]), "audited_full_span")

    def test_single_read_is_not_audited_full(self):
        self.assertEqual(self.cls([full("r1")]), "single_unit_full_span")

    def test_reverse_strand_counts(self):
        self.assertEqual(self.cls([full("r1", strand="-"), full("r2", off=7, strand="-")]), "audited_full_span")

    def test_x_deleted(self):
        deleted = lambda q: paf(q, 4000, 7500, "1000M1200D1300M")
        self.assertEqual(self.cls([deleted("d1"), deleted("d2")]), "conflict_x_absent")
        self.assertEqual(self.cls([full("r1"), full("r2", off=7), deleted("d1"), paf("d2", 3990, 7500, "1010M1200D1300M")]),
                         "audited_full_span_with_conflict")

    def test_secondary_alignment_lacking_x_is_not_conflict(self):
        # Two adjacent copies: a read's runner-up placement can skip X; that is not an absent X.
        secondary_skip = lambda q: paf(q, 4000, 7500, "1000M1200D1300M", tp="S", mapq=0)
        self.assertEqual(self.cls([full("r1"), full("r2", off=7), secondary_skip("s1"), secondary_skip("s2")]),
                         "audited_full_span")
        self.assertEqual(self.cls([secondary_skip("s1")]), "ambiguous_competing")

    def test_accumulated_small_indels_fail(self):
        cg = "1000M" + "17M3D" * 60 + "1300M"  # 180 bp of deletions inside a 1200 bp X
        self.assertEqual(self.cls([paf("r1", 4000, 7500, cg), paf("r2", 4000, 7500, cg)]), "full_span_cigar_fail")

    def test_missing_anchor_is_one_side(self):
        self.assertEqual(self.cls([paf("r1", 4950, 7500, "2550M"), paf("r2", 5100, 7500, "2400M")]), "one_side_only")

    def test_cigar_N_fails(self):
        self.assertEqual(self.cls([paf("r1", 4000, 7500, "1500M50N1950M"), paf("r2", 4000, 7500, "1500M50N1950M")]),
                         "full_span_cigar_fail")

    def test_equal_competing_position(self):
        lines = [full("r1"), paf("r1", 10000, 13500, "3500M", tp="S", t="c2", AS=7000),
                 full("r2"), paf("r2", 10000, 13500, "3500M", tp="S", t="c2", AS=6990)]
        self.assertEqual(self.cls(lines), "ambiguous_competing")

    def test_secondary_at_candidate(self):
        self.assertEqual(self.cls([full("r1", tp="S"), full("r2", tp="S")]), "ambiguous_competing")

    def test_equal_competing_adjacent_copy_on_same_contig(self):
        # Shift by 1300 bp: distinct position, inside the old L+200 exemption.
        lines = []
        for q, off in (("r1", 0), ("r2", 7)):
            lines += [full(q, off=off),
                      paf(q, 5300-off, 8800, f"{3500+off}M", tp="S", AS=7000, mapq=0)]
        audit = self.run_case(lines)["m1"]
        self.assertEqual(audit["audited_class"], "ambiguous_competing")
        self.assertEqual(audit["units_full_support"], "0")

    def test_identical_alignment_path_is_not_a_competing_position(self):
        lines = [full("r1"), full("r1", tp="S", mapq=0),
                 full("r2", off=7), full("r2", off=7, tp="S", mapq=0)]
        self.assertEqual(self.cls(lines), "audited_full_span")

    def test_secondary_cap_needs_recheck(self):
        lines = [full("r1"), full("r2", off=7)]
        for q in ("r1", "r2"):
            lines += [paf(q, 100 * i, 100 * i + 900, "900M", tp="S", t="c2", AS=1000) for i in range(50)]
        self.assertEqual(self.cls(lines), "ambiguous_competing")
        resolved = [full("r1"), full("r2", off=7)] + [paf(q, 100, 1000, "900M", tp="S", t="c2", AS=1000) for q in ("r1", "r2")]
        self.assertEqual(self.cls(lines, recheck=resolved), "audited_full_span")

    def test_duplicate_derived_read_is_one_unit(self):
        audit = self.run_case([full("r1"), full("r1_copy")])["m1"]
        self.assertEqual(audit["audited_class"], "single_unit_full_span")
        self.assertEqual(audit["duplicate_reads"], "1")

    def test_two_one_sided_reads_do_not_make_full_span(self):
        self.assertEqual(self.cls([paf("r1", 4000, 5500, "1500M"), paf("r2", 5800, 7500, "1700M")]), "one_side_only")

    def test_contig_edge_is_edge_limited(self):
        audit = self.run_case([paf("r1", 0, 3000, "3000M"), paf("r2", 0, 3000, "3000M")], members={"m1": EDGE})
        self.assertEqual(audit["m1"]["audited_class"], "edge_limited")

    def test_internal_only(self):
        self.assertEqual(self.cls([paf("r1", 5200, 6000, "800M")]), "internal_only")

    def test_checker_rejects_tampered_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "ref.fa").write_text(">c1\n" + "A" * 20000 + "\n")
            (tmp / "members.json").write_text(json.dumps([dict(member="m1", family="F", method="RS", accession="c1",
                candidate_start0=X[0], candidate_end0=X[1], boundary_source="provisional_homology_median")]))
            (tmp / "reads.paf").write_text(full("r1") + "\n")
            subprocess.run([sys.executable, str(SCRIPT), "classify", "--members", str(tmp / "members.json"),
                            "--reference", str(tmp / "ref.fa"), "--paf", str(tmp / "reads.paf"), "--sample", "t",
                            "--out", str(tmp / "out")], check=True, capture_output=True)
            audit = tmp / "out/member_audit.tsv"
            audit.write_text(audit.read_text().replace("single_unit_full_span", "audited_full_span"))
            check = subprocess.run([sys.executable, str(SCRIPT.with_name("check_audited_evidence.py")),
                                    "--audit-dir", str(tmp / "out"), "--members", str(tmp / "members.json"),
                                    "--reference", str(tmp / "ref.fa"), "--out", str(tmp / "check.json")],
                                   capture_output=True, text=True)
            self.assertNotEqual(check.returncode, 0)

    def test_legitimate_zero_and_malformed_input(self):
        self.assertEqual(self.cls([]), "no_coverage")
        self.run_case(["not\ta\tpaf"], expect_rc=1)


if __name__ == "__main__":
    unittest.main()
