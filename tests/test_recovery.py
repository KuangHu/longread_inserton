import json
import random
from dataclasses import replace

import pysam
import pytest

from site_recovery.core import Candidate, Config, aligned_slice, candidates_from_vcf, recover, write_outputs
from site_recovery.cli import execute, parser, validate_reference


def dna(n, seed):
    r = random.Random(seed)
    return "".join(r.choices("ACGT", k=n))


def fixture(tmp_path, kind="INS", other=3, bad=None, repeat=False):
    left, right, insert = dna(2000, 1), dna(2000, 2), dna(400, 3)
    insert = ("A" if right[0] != "A" else "C") + insert[1:-1] + ("G" if left[-1] != "G" else "T")
    if repeat:
        insert = insert[:-1]+left[-1]
    reference = left + (insert if kind == "DEL" else "") + right
    fa = tmp_path/"ref.fa"
    fa.write_text(">ctg\n"+reference+"\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "ctg", "LN": len(reference)}]}
    path = tmp_path/"reads.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i in range(4+other):
            event = i < 4
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"r{i}"
            read.flag = 16 if i % 2 else 0
            read.reference_id, read.reference_start, read.mapping_quality = 0, 1500, 60
            if event:
                read.cigartuples = [(0, 500), (1 if kind == "INS" else 2, 400), (0, 500)]
                read.query_sequence = left[-500:]+(insert if kind == "INS" else "")+right[:500]
                if bad == "low_mapq":
                    read.mapping_quality = 0
                elif bad == "softclip":
                    read.cigartuples = [(4, 500), (1 if kind == "INS" else 2, 400), (0, 500)]
            else:
                length = 1000+(400 if kind == "DEL" else 0)
                read.cigartuples = [(0, length)]
                read.query_sequence = reference[1500:1500+length]
                if bad == "skip":
                    read.cigartuples = [(0, 450), (3, 100), (0, length-550)]
                    read.query_sequence = reference[1500:1950]+reference[2050:1500+length]
            bam.write(read)
    pysam.index(str(path))
    c = Candidate("call1", "ctg", 1999, 2400 if kind == "DEL" else 2000, kind, 400)
    return fa, path, c, insert, left[-1000:]+right[:1000]


def run(fa, path, candidate, cfg=Config()):
    with pysam.FastaFile(str(fa)) as ref, pysam.AlignmentFile(str(path), "rb") as bam:
        return recover(candidate, bam, ref, cfg, "testhash", "sample")


@pytest.mark.parametrize("platform", ["hifi", "ont"])
@pytest.mark.parametrize("kind", ["INS", "DEL"])
def test_exact_recovery_both_routes(tmp_path, platform, kind):
    fa, path, c, insertion, empty = fixture(tmp_path, kind)
    rec = run(fa, path, c, Config(platform=platform))
    assert rec["status"] == "read_supported_pair"
    assert rec["sequences"]["insertion"] == insertion
    assert rec["sequences"]["empty_site"] == empty
    assert rec["empty_site_insertion_pos0"] == 1000
    assert rec["reference"]["start0"] == 2000
    assert rec["reference"]["end0"] == (2400 if kind == "DEL" else 2000)
    assert rec["support"]["n_occupied"] == (3 if kind == "DEL" else 4)
    assert rec["support"]["n_empty"] == (4 if kind == "DEL" else 3)


@pytest.mark.parametrize("kind", ["INS", "DEL"])
def test_no_other_allele_not_accepted(tmp_path, kind):
    fa, path, c, *_ = fixture(tmp_path, kind, other=0)
    assert run(fa, path, c)["status"] == "reference_backed_only"


@pytest.mark.parametrize("bad", ["low_mapq", "softclip"])
def test_invalid_anchors_not_recovered(tmp_path, bad):
    fa, path, c, *_ = fixture(tmp_path, bad=bad)
    assert run(fa, path, c)["status"] == "unresolved"


def test_reference_skip_not_empty_support(tmp_path):
    fa, path, c, *_ = fixture(tmp_path, bad="skip")
    rec = run(fa, path, c)
    assert rec["status"] == "reference_backed_only"
    assert rec["support"]["n_empty"] == 0


def test_repeat_reports_ambiguous(tmp_path):
    fa, path, c, *_ = fixture(tmp_path, repeat=True)
    rec = run(fa, path, c)
    assert rec["status"] == "ambiguous"
    assert rec["breakpoint_equivalence"]["start_min0"] < 2000


@pytest.mark.parametrize("change,reason", [({"svtype": "BND"}, "unsupported_svtype"),
                                         ({"filters": ("LowQual",)}, "vcf_filter_not_pass"),
                                         ({"svlen": 0}, "unknown_or_out_of_scope_length")])
def test_bad_candidates(tmp_path, change, reason):
    fa, path, c, *_ = fixture(tmp_path)
    assert run(fa, path, replace(c, **change))["reasons"] == [reason]


def test_cap_and_partial(tmp_path):
    fa, path, c, *_ = fixture(tmp_path)
    assert run(fa, path, c, Config(max_reads=2))["reasons"] == ["read_cap_exceeded_no_downsampled_call"]
    assert run(fa, path, c, Config(flank=3000))["status"] == "partial"


def test_dedup_outputs(tmp_path):
    fa, path, c, *_ = fixture(tmp_path)
    rec = run(fa, path, c)
    duplicate = run(fa, path, replace(c, source_id="call2"))
    out = tmp_path/"out"
    assert write_outputs([rec, duplicate], out) == {"read_supported_pair": 1}
    assert json.loads((out/"records.jsonl").read_text())["source_call_ids"] == ["call1", "call2"]
    assert (out/"accepted/insertion.fasta").read_text().count(">") == 1


def test_symbolic_vcf_recovers_cigar_not_host(tmp_path):
    fa, path, _, insertion, _ = fixture(tmp_path)
    vcf = tmp_path/"calls.vcf"
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=ctg,length=4000>\n'
                   '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
                   '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                   'ctg\t2000\tsymbolic\tA\t<INS>\t60\tPASS\tSVTYPE=INS;SVLEN=400\n'
                   'ctg\t2000\tbnd\tA\tA]ctg:3000]\t60\tPASS\tSVTYPE=BND\n')
    candidates = list(candidates_from_vcf(vcf))
    assert run(fa, path, candidates[0])["sequences"]["insertion"] == insertion
    assert run(fa, path, candidates[1])["status"] == "unresolved"
    args = parser().parse_args(["recover", "--reference", str(fa), "--bam", str(path),
                               "--vcf", str(vcf), "--platform", "ccs", "--sample", "x",
                               "--out", str(tmp_path/"result")])
    assert execute(args) == {"read_supported_pair": 1, "unresolved": 1}
    with pytest.raises(ValueError, match="new or empty"):
        execute(args)


def test_header_mismatch(tmp_path):
    fa, path, *_ = fixture(tmp_path)
    other = tmp_path/"other.fa"
    other.write_text(">ctg\nACGT\n")
    pysam.faidx(str(other))
    with pysam.FastaFile(str(other)) as ref, pysam.AlignmentFile(str(path), "rb") as bam:
        with pytest.raises(ValueError, match="exactly match"):
            validate_reference(bam, ref)


@pytest.mark.parametrize("cigar,sequence,expected", [
    ([(0, 100), (1, 2), (0, 100)], "A"*100+"CC"+"A"*100, "A"*50+"CC"+"A"*50),
    ([(0, 100), (2, 2), (0, 100)], "A"*200, "A"*98),
    ([(0, 100), (2, 50), (0, 100)], "A"*200, None),
    ([(0, 100), (3, 2), (0, 100)], "A"*200, None),
])
def test_small_errors_tolerated_but_sv_and_skips_rejected(cigar, sequence, expected):
    read = pysam.AlignedSegment()
    read.reference_start = 0
    read.query_sequence = sequence
    read.cigartuples = cigar
    assert aligned_slice(read, 50, 150) == expected


def test_distant_locus_not_recruited_by_insert_sequence(tmp_path):
    fa, path, c, *_ = fixture(tmp_path)
    # Same reference and insert length, different locus: supporting insert
    # sequence alone must not create a locus call.
    rec = run(fa, path, replace(c, start0=3000, end0=3001))
    assert rec["status"] == "unresolved"
    assert not rec["sequences"]


def test_supplementary_only_not_promoted(tmp_path):
    fa, path, c, *_ = fixture(tmp_path)
    new = tmp_path/"supp.bam"
    with pysam.AlignmentFile(str(path), "rb") as source, pysam.AlignmentFile(str(new), "wb", header=source.header) as out:
        for read in source:
            if read.query_name in {"r0", "r1", "r2", "r3"}:
                read.flag |= 2048
            out.write(read)
    pysam.index(str(new))
    assert run(fa, new, c)["status"] == "unresolved"
