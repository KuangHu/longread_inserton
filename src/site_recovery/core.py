"""Conservative, locus-anchored recovery from full CIGAR INS/DEL evidence.

VCF coordinates only nominate a search region. All reported breakpoints are
derived from CIGAR operations, on the input assembly's forward strand.
Split-only and multi-allelic/ambiguous cases are retained as unresolved.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
import hashlib
import json
from pathlib import Path
from typing import Iterator

import edlib
import pysam
import spoa

MATCH = {0, 7, 8}
REF_CONSUME = {0, 2, 3, 7, 8}
QUERY_CONSUME = {0, 1, 4, 7, 8}
DNA = set("ACGT")


@dataclass(frozen=True)
class Config:
    platform: str = "hifi"
    flank: int = 1000
    anchor: int = 100
    min_mapq: int = 20
    min_support: int = 3
    min_other_support: int = 2
    min_size: int = 50
    max_size: int = 100000
    search_slop: int = 100
    length_tolerance: float = 0.25
    max_sequences: int = 100
    max_reads: int = 10000
    repeat_limit: int = 100

    @property
    def identity(self) -> float:
        return 0.95 if self.platform == "hifi" else 0.90

    def validate(self):
        if self.platform not in {"hifi", "ont"}:
            raise ValueError("platform must be hifi or ont")
        for key in ("flank", "anchor", "min_support", "min_other_support", "min_size",
                    "max_size", "max_sequences", "max_reads", "repeat_limit"):
            if getattr(self, key) < 1:
                raise ValueError(f"{key} must be positive")
        if self.flank < self.anchor or self.max_size < self.min_size:
            raise ValueError("flank must be >= anchor, and max_size >= min_size")
        if self.min_mapq < 1 or self.search_slop < 0:
            raise ValueError("min_mapq must be positive; search_slop must be nonnegative")
        if not 0 <= self.length_tolerance < 1:
            raise ValueError("length_tolerance must lie in [0,1)")


@dataclass(frozen=True)
class Candidate:
    source_id: str
    chrom: str
    start0: int
    end0: int
    svtype: str
    svlen: int
    filters: tuple[str, ...] = ()


@dataclass(frozen=True)
class Gap:
    read_name: str
    start0: int
    end0: int
    qstart: int
    qend: int
    index: int
    op: int
    sequence: str


def digest(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return 1 - edlib.align(a, b, mode="NW", task="distance")["editDistance"] / max(len(a), len(b))


def candidates_from_vcf(path: Path) -> Iterator[Candidate]:
    with pysam.VariantFile(str(path)) as vcf:
        for i, rec in enumerate(vcf):
            kind = str(rec.info.get("SVTYPE", "UNKNOWN"))
            raw = rec.info.get("SVLEN", 0)
            if isinstance(raw, (tuple, list)):
                raw = raw[0] if len(raw) == 1 else 0
            size = abs(int(raw or 0))
            # Never parse BND ALT syntax as sequence. Never infer an unbounded
            # symbolic INS length from a host contig.
            if not size and kind == "DEL":
                size = max(0, rec.stop - rec.start - 1)
            if not size and kind == "INS" and rec.alts and len(rec.alts) == 1:
                alt = rec.alts[0].upper()
                if set(alt) <= DNA:
                    size = max(0, len(alt) - len(rec.ref))
            if rec.alts and len(rec.alts) > 1:
                kind = "MULTIALLELIC"
            yield Candidate(rec.id or f"row_{i+1}", rec.chrom, rec.start,
                            rec.stop, kind, size, tuple(rec.filter.keys()))


def usable_read(read, cfg: Config) -> bool:
    return (not read.is_unmapped and not read.is_secondary and
            not read.is_supplementary and not read.is_duplicate and
            not read.is_qcfail and read.mapping_quality >= cfg.min_mapq and
            read.query_sequence is not None)


def cigar_gaps(read) -> Iterator[Gap]:
    rpos, qpos = read.reference_start, 0
    seq = read.query_sequence or ""
    for i, (op, n) in enumerate(read.cigartuples or []):
        if op == 1:
            yield Gap(read.query_name, rpos, rpos, qpos, qpos+n, i, op, seq[qpos:qpos+n].upper())
        elif op == 2:
            yield Gap(read.query_name, rpos, rpos+n, qpos, qpos, i, op, "")
        if op in REF_CONSUME:
            rpos += n
        if op in QUERY_CONSUME:
            qpos += n


def anchored_gap(read, gap: Gap, ref: pysam.FastaFile, cfg: Config) -> bool:
    """Both anchors must be real query bases, not D/N spans or soft clipping."""
    a = cfg.anchor
    if gap.start0 < a or gap.end0+a > ref.get_reference_length(read.reference_name):
        return False
    query_left = aligned_slice(read, gap.start0-a, gap.start0)
    query_right = aligned_slice(read, gap.end0, gap.end0+a)
    if query_left is None or query_right is None:
        return False
    left = ref.fetch(read.reference_name, gap.start0-a, gap.start0).upper()
    right = ref.fetch(read.reference_name, gap.end0, gap.end0+a).upper()
    return (set(left+right) <= DNA and
            identity(query_left, left) >= cfg.identity and
            identity(query_right, right) >= cfg.identity)


def aligned_slice(read, start0: int, end0: int, max_small_indel: int = 10) -> str | None:
    """Extract an anchored query interval, allowing bounded sequencing indels.

    Both boundary bases must align to real query bases. N and internal I/D
    longer than max_small_indel invalidate support. Identity is checked by
    the caller. An insertion exactly at an interval boundary is excluded.
    """
    rpos, qpos = read.reference_start, 0
    qstart = qend = None
    for op, n in read.cigartuples or []:
        if op in MATCH:
            if rpos <= start0 < rpos+n:
                qstart = qpos+start0-rpos
            if rpos < end0 <= rpos+n:
                qend = qpos+end0-rpos
        elif op == 1 and start0 < rpos < end0 and n > max_small_indel:
            return None
        elif op in {2, 3} and rpos < end0 and rpos+n > start0:
            if op == 3 or n > max_small_indel:
                return None
        if op in REF_CONSUME:
            rpos += n
        if op in QUERY_CONSUME:
            qpos += n
    if qstart is None or qend is None or qend <= qstart:
        return None
    return read.query_sequence[qstart:qend].upper()


def repeat_equivalence(ref: pysam.FastaFile, chrom: str, start: int, end: int,
                       insertion: str, limit: int) -> dict:
    """Report exact repeat-induced breakpoint ambiguity; never silently fix it.

    Positions are alternative *starts*, in 0-based interbase coordinates.
    A deleted interval has the same length at every equivalent start.
    """
    n = ref.get_reference_length(chrom)
    left, right = start, start
    ins_left = ins_right = insertion
    width = end-start
    for _ in range(limit):
        if left == 0:
            break
        before = ref.fetch(chrom, left-1, left).upper()
        tail = (ref.fetch(chrom, left+width-1, left+width).upper()
                if width else ins_left[-1:])
        if before not in DNA or before != tail:
            break
        left -= 1
        if not width:
            ins_left = ins_left[-1]+ins_left[:-1]
    for _ in range(limit):
        if right+width >= n:
            break
        after = ref.fetch(chrom, right+width, right+width+1).upper()
        head = ref.fetch(chrom, right, right+1).upper() if width else ins_right[:1]
        if after not in DNA or after != head:
            break
        right += 1
        if not width:
            ins_right = ins_right[1:]+ins_right[0]
    return {"start_min0": left, "start_max0": right,
            "range_may_extend": start-left == limit or right-start == limit,
            "ambiguous": left != right}


def choose_gap_group(pairs, cfg: Config):
    """Cluster exact deletion intervals / insertion points; preserve alternatives."""
    groups = defaultdict(list)
    for read, gap in pairs:
        groups[(gap.start0, gap.end0)].append((read, gap))
    ranked = sorted(groups.items(), key=lambda item: (-len(item[1]), item[0]))
    if not ranked:
        return [], []
    return ranked[0][1], [{"start0": key[0], "end0": key[1], "reads": len(value)}
                         for key, value in ranked]


def recover(candidate: Candidate, bam: pysam.AlignmentFile, ref: pysam.FastaFile,
            cfg: Config, reference_sha256: str, sample: str) -> dict:
    result = {"source_call": asdict(candidate), "sample": sample,
              "status": "unresolved", "reasons": [], "platform": cfg.platform,
              "coordinate_system": "0-based; intervals half-open; positions interbase",
              "orientation": "assembly_forward", "sequences": {}}
    kind = candidate.svtype
    if kind not in {"INS", "DEL"}:
        result["reasons"] = ["unsupported_svtype"]
        return result
    if set(candidate.filters)-{"PASS", "."}:
        result["reasons"] = ["vcf_filter_not_pass"]
        return result
    if not cfg.min_size <= candidate.svlen <= cfg.max_size:
        result["reasons"] = ["unknown_or_out_of_scope_length"]
        return result
    if candidate.chrom not in ref.references or candidate.chrom not in bam.references:
        result["reasons"] = ["missing_contig"]
        return result
    n = ref.get_reference_length(candidate.chrom)
    lower = max(0, candidate.start0-cfg.search_slop-cfg.anchor)
    upper = min(n, candidate.start0+cfg.search_slop+cfg.anchor+1+
                (candidate.svlen if kind == "DEL" else 0))
    if lower >= upper:
        result["reasons"] = ["candidate_out_of_bounds"]
        return result
    reads = []
    counts = Counter()
    name_counts = Counter()
    for read in bam.fetch(candidate.chrom, lower, upper):
        counts["fetched"] += 1
        if counts["fetched"] > cfg.max_reads:
            result["reasons"] = ["read_cap_exceeded_no_downsampled_call"]
            return result
        if not usable_read(read, cfg):
            counts["excluded_alignment"] += 1
            continue
        reads.append(read)
        name_counts[read.query_name] += 1
    # Multiple primary records bearing one name are ambiguous molecules.
    reads = [r for r in reads if name_counts[r.query_name] == 1]
    counts["usable_unique_primary"] = len(reads)
    expected_op = 1 if kind == "INS" else 2
    length_slop = max(5, round(candidate.svlen*cfg.length_tolerance))
    pairs = []
    for read in reads:
        hits = []
        for gap in cigar_gaps(read):
            length = gap.qend-gap.qstart if gap.op == 1 else gap.end0-gap.start0
            if (gap.op == expected_op and abs(gap.start0-candidate.start0) <= cfg.search_slop
                    and abs(length-candidate.svlen) <= length_slop
                    and anchored_gap(read, gap, ref, cfg)):
                hits.append(gap)
        if len(hits) == 1:
            pairs.append((read, hits[0]))
    group, alternatives = choose_gap_group(pairs, cfg)
    result["read_counts"] = dict(counts)
    result["breakpoint_groups"] = alternatives
    if not group:
        result["reasons"] = ["no_full_cigar_event_with_two_anchors"]
        return result
    start, end = group[0][1].start0, group[0][1].end0
    result["reference"] = {"contig": candidate.chrom, "sha256": reference_sha256,
                            "start0": start, "end0": end,
                            "occupied_in_reference": kind == "DEL"}
    if len(group) < cfg.min_support:
        result["reasons"] = ["insufficient_same_breakpoint_support"]
        return result
    competing = any(g["reads"] >= cfg.min_support for g in alternatives[1:])
    if competing:
        result["reasons"] = ["multiple_supported_breakpoints"]
        return result
    chosen = sorted(group, key=lambda p: (-p[0].mapping_quality, p[0].query_name))[:cfg.max_sequences]
    if kind == "INS":
        sequences = [gap.sequence for _, gap in chosen]
        if any(not s or set(s)-DNA for s in sequences):
            result["reasons"] = ["ambiguous_insertion_bases"]
            return result
        insertion, _ = spoa.poa(sequences, algorithm=1)
        insertion = insertion.upper()
        # Do not average distinct insert alleles into a chimeric consensus.
        agree = [identity(gap.sequence, insertion) >= cfg.identity for _, gap in group]
        if not all(agree):
            result["reasons"] = ["discordant_insertion_sequences"]
            return result
        sequence_source = "POA_of_locus_anchored_CIGAR_insertions"
    else:
        insertion = ref.fetch(candidate.chrom, start, end).upper()
        sequence_source = "assembly_interval_with_deletion_read_support"
    if not cfg.min_size <= len(insertion) <= cfg.max_size or set(insertion)-DNA:
        result["reasons"] = ["invalid_recovered_insertion"]
        return result
    left_start, right_end = max(0, start-cfg.flank), min(n, end+cfg.flank)
    upstream = ref.fetch(candidate.chrom, left_start, start).upper()
    downstream = ref.fetch(candidate.chrom, end, right_end).upper()
    if set(upstream+downstream)-DNA:
        result["reasons"] = ["ambiguous_flank_bases"]
        return result
    # Reference-state support must bridge the entire interval with actual
    # sequence identity, tolerating small sequencing indels but no SV gaps.
    ref_start, ref_end = start-cfg.anchor, end+cfg.anchor
    ref_allele = ref.fetch(candidate.chrom, ref_start, ref_end).upper()
    event_names = {g.read_name for _, g in group}
    other_names = []
    for read in reads:
        if read.query_name in event_names:
            continue
        seq = aligned_slice(read, ref_start, ref_end)
        if seq is not None and identity(seq, ref_allele) >= cfg.identity:
            other_names.append(read.query_name)
    occupied_names = sorted(event_names) if kind == "INS" else sorted(other_names)
    empty_names = sorted(other_names) if kind == "INS" else sorted(event_names)
    result["support"] = {"occupied_reads": occupied_names, "empty_reads": empty_names,
                         "n_occupied": len(occupied_names), "n_empty": len(empty_names),
                         "consensus_read_names": [g.read_name for _, g in chosen] if kind == "INS" else [],
                         "evidence_scope": "local; does not phase a whole strain genome",
                         "anchor_bp": cfg.anchor, "min_mapq": cfg.min_mapq,
                         "min_sequence_identity": cfg.identity}
    result["sequences"] = {"insertion": insertion, "upstream": upstream,
                           "downstream": downstream, "empty_site": upstream+downstream}
    result["insertion_sequence_source"] = sequence_source
    result["empty_site_sequence_source"] = "assembly_flanks; no strain-specific SNP polishing"
    result["empty_site_insertion_pos0"] = len(upstream)
    result["empty_site_flanks"] = {"upstream_reference_start0": left_start,
                                  "upstream_reference_end0": start,
                                  "downstream_reference_start0": end,
                                  "downstream_reference_end0": right_end}
    result["full_requested_flanks"] = len(upstream) == len(downstream) == cfg.flank
    equiv = repeat_equivalence(ref, candidate.chrom, start, end, insertion, cfg.repeat_limit)
    result["breakpoint_equivalence"] = equiv
    if equiv["ambiguous"]:
        result["status"] = "ambiguous"
        result["reasons"] = ["repeat_equivalent_breakpoints"]
    elif len(other_names) < cfg.min_other_support:
        result["status"] = "reference_backed_only"
        result["reasons"] = ["other_state_not_sufficiently_read_supported"]
    elif not result["full_requested_flanks"]:
        result["status"] = "partial"
        result["reasons"] = ["requested_flanks_truncated_at_contig_edge"]
    else:
        result["status"] = "read_supported_pair"
    # No random UUIDs; distinguish different references and insert alleles.
    event_key = [reference_sha256, candidate.chrom, start, end, digest(insertion)]
    result["event_id"] = "isr_"+digest(json.dumps(event_key, separators=(",", ":")))[:20]
    return result


def write_outputs(records: list[dict], out: Path):
    import csv
    out.mkdir(parents=True, exist_ok=True)
    # Deduplicate repeated VCF observations; keep all source IDs in JSON.
    unique = {}
    unresolved = []
    for rec in records:
        if "event_id" not in rec:
            unresolved.append(rec)
        elif rec["event_id"] not in unique:
            rec["source_call_ids"] = [rec["source_call"]["source_id"]]
            unique[rec["event_id"]] = rec
        else:
            unique[rec["event_id"]]["source_call_ids"].append(rec["source_call"]["source_id"])
    records = list(unique.values())+unresolved
    with (out/"records.jsonl").open("w") as fh:
        for rec in records:
            fh.write(json.dumps(rec, ensure_ascii=False)+"\n")
    fields = ["event_id", "sample", "status", "reasons", "source_call_ids", "svtype",
              "contig", "reference_start0", "reference_end0", "insertion_length",
              "upstream_length", "downstream_length", "empty_site_insertion_pos0",
              "occupied_read_count", "empty_read_count"]
    with (out/"events.tsv").open("w") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rec in records:
            ref = rec.get("reference", {})
            seq = rec["sequences"]
            sup = rec.get("support", {})
            writer.writerow(dict(event_id=rec.get("event_id", ""), sample=rec["sample"],
                status=rec["status"], reasons=";".join(rec["reasons"]),
                source_call_ids=";".join(rec.get("source_call_ids", [rec["source_call"]["source_id"]])),
                svtype=rec["source_call"]["svtype"], contig=rec["source_call"]["chrom"],
                reference_start0=ref.get("start0", ""), reference_end0=ref.get("end0", ""),
                insertion_length=len(seq.get("insertion", "")), upstream_length=len(seq.get("upstream", "")),
                downstream_length=len(seq.get("downstream", "")),
                empty_site_insertion_pos0=rec.get("empty_site_insertion_pos0", ""),
                occupied_read_count=sup.get("n_occupied", ""), empty_read_count=sup.get("n_empty", "")))
    for folder, selected in (("accepted", [r for r in records if r["status"] == "read_supported_pair"]),
                             ("candidates", [r for r in records if r["sequences"] and r["status"] != "read_supported_pair"])):
        target = out/folder
        target.mkdir(exist_ok=True)
        for key in ("insertion", "empty_site", "upstream", "downstream"):
            with (target/f"{key}.fasta").open("w") as fh:
                for rec in selected:
                    seq = rec["sequences"][key]
                    fh.write(f'>{rec["event_id"]} status={rec["status"]} insertion_pos0={rec["empty_site_insertion_pos0"]}\n')
                    for i in range(0, len(seq), 80):
                        fh.write(seq[i:i+80]+"\n")
    summary = dict(Counter(r["status"] for r in records))
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    return summary
