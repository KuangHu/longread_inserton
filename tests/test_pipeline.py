import json
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from multisite_insertion.boundaries import core_link_count, infer
from multisite_insertion.io import load_manifest, read_json, write_json, fasta_dict
from multisite_insertion.native import overlap_parts, collect
from multisite_insertion.pipeline import prepare_contexts, run, select_families


class BoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out = Path(self.temp.name)
        for folder in ["inputs", "results", "sequences"]:
            (self.out / folder).mkdir()
        self.genomes = {"sample": {"c": SeqRecord(Seq("ACGT" * 1250), id="c")}}
        self.family = dict(
            query="sample_PILER_0.0",
            organism="sample",
            method="PILER",
            raw_occurrences=3,
            best_deposited_tnp_hit={"title": "transposase"},
            native_records_overlapping_tnp_CDS=1,
            occurrences=[
                dict(
                    contig="c",
                    accession="c",
                    start0=a,
                    end0=a + 100,
                    strand=strand,
                    in_bounds=True,
                )
                for a, strand in [(1000, "+"), (2000, "-"), (3000, "+")]
            ],
        )
        self.selected = {self.family["query"]: self.family}

    def contexts(self, hits=True):
        meta, invalid = prepare_contexts(self.selected, self.genomes, self.out)
        lines = []
        if hits:
            for q, m in meta.items():
                for s, n in meta.items():
                    if q == s:
                        continue
                    lines.append(
                        "\t".join(
                            map(
                                str,
                                [
                                    q,
                                    s,
                                    100,
                                    100,
                                    m["core_start0"] + 1,
                                    m["core_end0"],
                                    n["core_start0"] + 1,
                                    n["core_end0"],
                                    "1e-40",
                                    200,
                                    m["extract_end0"] - m["extract_start0"],
                                    n["extract_end0"] - n["extract_start0"],
                                    "100",
                                ],
                            )
                        )
                    )
        (self.out / "results/member_contexts_btop.blastn.tsv").write_text(
            "\n".join(lines) + ("\n" if lines else "")
        )
        return meta, invalid

    def test_reverse_orientation_round_trip(self):
        meta, _ = self.contexts()
        rows = infer(self.selected, self.genomes, {"sample": []}, self.out, meta)
        r = rows[1]
        self.assertEqual((r["candidate_start0"], r["candidate_end0"]), (2000, 2100))
        self.assertEqual(r["upstream_genomic_start0"], 2100)
        exported = fasta_dict(self.out / "sequences/triplet.fasta")
        self.assertEqual(
            str(exported[r["member"]].seq),
            str(self.genomes["sample"]["c"].seq[1500:2600].reverse_complement()),
        )

    def test_adjacent_annotation_two_bp_is_not_adopted(self):
        meta, _ = self.contexts()
        feature = dict(
            id="adjacent",
            accession="c",
            type="mobile_element",
            start0=1098,
            end0=1978,
            parts=[dict(start0=1098, end0=1978)],
            qualifiers={"mobile_element_type": ["insertion sequence:test"]},
        )
        row = infer(self.selected, self.genomes, {"sample": [feature]}, self.out, meta)[
            0
        ]
        self.assertEqual(
            row["overlapping_mobile_annotations"][0]["overlap_native_bp"], 2
        )
        self.assertEqual(row["candidate_length"], 100)
        self.assertEqual(row["boundary_source"], "provisional_homology_median")

    def test_annotation_and_inferred_ranges_remain_separate(self):
        meta, _ = self.contexts()
        feature = dict(
            id="element",
            accession="c",
            type="mobile_element",
            start0=1001,
            end0=1099,
            parts=[dict(start0=1001, end0=1099)],
            qualifiers={},
        )
        row = infer(self.selected, self.genomes, {"sample": [feature]}, self.out, meta)[
            0
        ]
        self.assertEqual((row["inferred_start0"], row["inferred_end0"]), (1000, 1100))
        self.assertEqual((row["candidate_start0"], row["candidate_end0"]), (1001, 1099))
        self.assertEqual(row["boundary_source"], "deposited_mobile_element")

    def test_compound_element_does_not_become_whole_chromosome(self):
        meta, _ = self.contexts()
        feature = dict(
            id="compound",
            accession="c",
            type="mobile_element",
            start0=0,
            end0=5000,
            parts=[dict(start0=0, end0=1090), dict(start0=4900, end0=5000)],
            qualifiers={},
        )
        row = infer(self.selected, self.genomes, {"sample": [feature]}, self.out, meta)[
            0
        ]
        self.assertEqual(row["candidate_length"], 100)
        self.assertEqual(row["boundary_source"], "provisional_homology_median")

    def test_invalid_record_retained_without_sequence(self):
        self.family["occurrences"][0].update(start0=-5, in_bounds=False)
        meta, invalid = self.contexts()
        self.assertEqual(len(meta), 2)
        self.assertEqual(len(invalid), 1)
        self.assertEqual(invalid[0]["native_record"]["start0"], -5)
        self.assertFalse(invalid[0]["sequence_exported"])

    def test_unresolved_and_empty_flank_are_not_discarded(self):
        self.family["occurrences"][0].update(start0=0, end0=100)
        meta, _ = self.contexts(hits=False)
        rows = infer(self.selected, self.genomes, {"sample": []}, self.out, meta)
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]["boundary_source"], "native_interval_only_unresolved")
        self.assertEqual(rows[0]["upstream_length"], 0)
        self.assertFalse(rows[0]["flanks_complete"])

    def test_window_edge_is_censored(self):
        meta, _ = self.contexts(hits=False)
        q, s = list(meta)[:2]
        # An HSP links both cores and reaches the entire query window start.
        m, n = meta[q], meta[s]
        offset = n["core_start0"] - m["core_start0"]
        self.assertGreaterEqual(offset, 0)
        length = m["core_end0"]
        line = "\t".join(
            map(
                str,
                [
                    q,
                    s,
                    100,
                    length,
                    1,
                    length,
                    offset + 1,
                    offset + length,
                    "1e-100",
                    1000,
                    m["extract_end0"] - m["extract_start0"],
                    n["extract_end0"] - n["extract_start0"],
                    length,
                ],
            )
        )
        (self.out / "results/member_contexts_btop.blastn.tsv").write_text(line + "\n")
        row = infer(self.selected, self.genomes, {"sample": []}, self.out, meta)[0]
        self.assertEqual(row["boundary_source"], "homology_extent_censored_at_window")

    def test_no_selected_families_produces_empty_files(self):
        self.selected = {}
        meta, _ = self.contexts(hits=False)
        self.assertEqual(infer({}, self.genomes, {"sample": []}, self.out, meta), [])
        self.assertEqual(
            (self.out / "sequences/candidate_insertion.fasta").read_text(), ""
        )

    def test_family_support_does_not_require_every_member(self):
        selected = select_families([self.family], "tnp-supported")
        self.assertEqual(len(selected[self.family["query"]]["occurrences"]), 3)
        self.family["best_deposited_tnp_hit"] = None
        self.assertEqual(select_families([self.family], "tnp-supported"), {})
        self.assertEqual(len(select_families([self.family], "all")), 1)


class DataTests(unittest.TestCase):
    def test_all_invalid_members_complete_without_blast_or_silent_loss(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "native").mkdir()
            (p / "g.fa").write_text(">c\n" + "ACGT" * 100 + "\n")
            (p / "native/families.fasta").write_text(">R=0\nACGT\n")
            (p / "native/ranges.tsv").write_text("c\t390\t420\tR=0\t+\n")
            write_json(
                p / "m.json",
                {
                    "schema_version": 1,
                    "samples": [
                        {"id": "s", "fasta": "g.fa", "native": {"RS": "native"}}
                    ],
                },
            )
            state = run(p / "m.json", p / "out", methods=["RS"], subset="all")
            self.assertEqual(state["status"], "complete")
            self.assertEqual(state["summary"]["selected_native_records"], 1)
            self.assertEqual(state["summary"]["invalid_records"], 1)
            self.assertEqual(state["summary"]["exported_records"], 0)
            self.assertEqual(len(read_json(p / "out/results/record_inventory.json")), 1)

    def test_empty_native_discovery_completes_without_external_tools(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "native").mkdir()
            (p / "g.fa").write_text(">c\nACGT\n")
            for f in ["families.fasta", "ranges.tsv"]:
                (p / "native" / f).write_text("")
            write_json(
                p / "m.json",
                {
                    "schema_version": 1,
                    "samples": [
                        {"id": "s", "fasta": "g.fa", "native": {"RS": "native"}}
                    ],
                },
            )
            state = run(p / "m.json", p / "out", methods=["RS"])
            self.assertEqual(state["status"], "complete")
            self.assertEqual(state["summary"]["native_families"], 0)
            self.assertEqual(
                (p / "out/sequences/candidate_insertion.fasta").read_text(), ""
            )

    def test_btop_gaps_consume_correct_axes(self):
        meta = dict(core_start0=0, core_end0=100)
        self.assertEqual(core_link_count(0, 0, "3AC2A-2-G4", meta, meta), (12, 13, 13))

    def test_overlapping_windows_outside_partner_core_do_not_support_boundaries(self):
        self.assertEqual(
            core_link_count(
                0,
                0,
                "1000",
                dict(core_start0=100, core_end0=200),
                dict(core_start0=700, core_end0=800),
            )[0],
            0,
        )

    def test_compound_origin_cds_overlap(self):
        parts = [dict(start0=950, end0=1000), dict(start0=0, end0=50)]
        self.assertEqual(overlap_parts(400, 600, parts), 0)
        self.assertEqual(overlap_parts(0, 30, parts), 30)

    def test_failure_is_not_zero_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "genome.fa").write_text(">c\n" + "ACGT" * 100 + "\n")
            write_json(
                p / "manifest.json",
                dict(schema_version=1, samples=[dict(id="one", fasta="genome.fa")]),
            )
            with self.assertRaises(ValueError):
                run(
                    p / "manifest.json",
                    p / "out",
                    methods=["PILER"],
                    bins={"pals": "/nonexistent/native_binary"},
                )
            self.assertEqual(read_json(p / "out/run.json")["status"], "failed")

    def test_duplicate_sample_id_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "g.fa").write_text(">c\nACGT\n")
            write_json(
                p / "m.json",
                dict(
                    schema_version=1,
                    samples=[dict(id="s", fasta="g.fa"), dict(id="s", fasta="g.fa")],
                ),
            )
            with self.assertRaises(ValueError):
                load_manifest(p / "m.json")

    def test_native_duplicate_and_out_of_bounds_occurrences_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            for d in ["inputs", "results", "native/s/RS"]:
                (p / d).mkdir(parents=True, exist_ok=True)
            (p / "inputs/g.fa").write_text(">c\n" + "ACGT" * 100 + "\n")
            (p / "native/s/RS/families.fasta").write_text(">R=0\nACGT\n")
            (p / "native/s/RS/ranges.tsv").write_text(
                "c\t5\t20\tR=0\t+\nc\t5\t20\tR=0\t+\nc\t390\t420\tR=0\t-\n"
            )
            rows = collect(
                [dict(id="s", fasta="inputs/g.fa")], p, ["RS"], {"s": []}, {"s": {}}
            )
            self.assertEqual(rows[0]["raw_occurrences"], 3)
            self.assertEqual(sum(o["in_bounds"] for o in rows[0]["occurrences"]), 2)
            self.assertEqual(rows[0]["occurrences"][2]["end0"], 420)


if __name__ == "__main__":
    unittest.main()
