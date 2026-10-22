"""One assembly set per sample; annotation and boundary inference follow native discovery."""

from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import __version__
from .annotation import annotate, prepare_annotations
from .boundaries import infer
from .io import (
    Runner,
    completed_blast,
    digest,
    fasta_dict,
    load_manifest,
    read_json,
    samples_in_run,
    stage_inputs,
    write_json,
)
from .native import collect, discover
from .summary import summarize

BLASTN_FORMAT = (
    "7 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"
)


def select_families(families, subset):
    if subset not in ("all", "tnp-supported"):
        raise ValueError("Unknown subset: " + subset)
    return {
        f["query"]: f
        for f in families
        if subset == "all"
        or (
            f["best_deposited_tnp_hit"] is not None
            and f["native_records_overlapping_tnp_CDS"] > 0
        )
    }


def prepare_contexts(selected, genomes, out):
    records, meta, invalid = [], {}, []
    for family in selected.values():
        for i, o in enumerate(family["occurrences"], 1):
            q = f'{family["query"]}_m{i}'
            if not o["in_bounds"]:
                invalid.append(
                    dict(
                        member=q,
                        family=family["query"],
                        organism=family["organism"],
                        method=family["method"],
                        boundary_source="invalid_native_coordinates",
                        sequence_exported=False,
                        native_record=o,
                        reason="Native coordinates retained verbatim; no clipping or invented sequence.",
                    )
                )
                continue
            seq = genomes[family["organism"]][o["contig"]].seq
            a, b = max(0, o["start0"] - 2000), min(len(seq), o["end0"] + 2000)
            s = seq[a:b]
            q0, q1 = o["start0"] - a, o["end0"] - a
            if o["strand"] == "-":
                s = s.reverse_complement()
                q0, q1 = len(s) - q1, len(s) - q0
            records.append(SeqRecord(s, id=q, description=""))
            meta[q] = dict(
                family=family["query"],
                method=family["method"],
                organism=family["organism"],
                contig=o["contig"],
                native_start0=o["start0"],
                native_end0=o["end0"],
                extract_start0=a,
                extract_end0=b,
                core_start0=q0,
                core_end0=q1,
            )
    SeqIO.write(records, Path(out) / "inputs/selected_member_contexts.fasta", "fasta")
    write_json(Path(out) / "inputs/member_contexts.json", meta)
    write_json(Path(out) / "results/invalid_native_records.json", invalid)
    return meta, invalid


def blast_pairs(runner, query, output, label, count, btop=False):
    if not count:
        Path(output).write_text("# No sequences to compare.\n")
        return
    fmt = BLASTN_FORMAT + (" btop" if btop else "")
    runner.run(
        label,
        "blastn",
        [
            "-query",
            query,
            "-subject",
            query,
            "-task",
            "blastn",
            "-word_size",
            "11",
            "-evalue",
            "1e-5",
            "-dust",
            "yes",
            "-soft_masking",
            "true",
            "-max_target_seqs",
            "100000",
            "-outfmt",
            fmt,
            "-out",
            output,
        ],
    )
    completed_blast(output, count)


def finish_boundaries(out, runner, families, genomes, features, subset):
    selected = select_families(families, subset)
    meta, invalid = prepare_contexts(selected, genomes, out)
    blast_pairs(
        runner,
        out / "inputs/selected_member_contexts.fasta",
        out / "results/member_contexts_btop.blastn.tsv",
        "member_contexts_blastn",
        len(meta),
        btop=True,
    )
    rows = infer(selected, genomes, features, out, meta)
    # Invalid records remain in a separate, explicitly counted result, never clipped.
    write_json(
        out / "results/record_inventory.json",
        [
            dict(
                member=r["member"],
                family=r["family"],
                sequence_exported=True,
                boundary_source=r["boundary_source"],
            )
            for r in rows
        ]
        + invalid,
    )
    flanks = out / "sequences/flanks_for_blast.fasta"
    flank_records = fasta_dict(flanks)
    count = sum(len(r) > 0 for r in flank_records.values())
    # BLAST accepts truncated flanks, but not zero-length FASTA entries.
    usable = out / "inputs/nonempty_flanks.fasta"
    SeqIO.write((r for r in flank_records.values() if len(r)), usable, "fasta")
    blast_pairs(
        runner,
        usable,
        out / "results/flanks_all_vs_all.blastn.tsv",
        "flanks_blastn",
        count,
    )
    summaries = summarize(out)
    write_json(
        out / "results/selected_family_inventory.json",
        [
            dict(
                family=f["query"],
                native_records=f["raw_occurrences"],
                exportable_records=sum(o["in_bounds"] for o in f["occurrences"]),
            )
            for f in selected.values()
        ],
    )
    totals = dict(
        package_version=__version__,
        subset=subset,
        native_families=len(families),
        native_records=sum(f["raw_occurrences"] for f in families),
        selected_families=len(selected),
        selected_native_records=sum(f["raw_occurrences"] for f in selected.values()),
        exported_records=len(rows),
        invalid_records=len(invalid),
        boundary_sources=dict(Counter(r["boundary_source"] for r in rows)),
        flanks_complete=sum(r["flanks_complete"] for r in rows),
        families_with_exportable_records=len(summaries),
    )
    write_json(out / "results/summary.json", totals)
    return totals


def run(
    manifest,
    out,
    methods=("PILER", "RS"),
    subset="tnp-supported",
    bins=None,
    threads=2,
    discovery_only=False,
):
    if threads < 1:
        raise ValueError("threads must be positive")
    samples = load_manifest(manifest)
    out = Path(out).resolve()
    samples = stage_inputs(samples, out)
    for d in ["results", "reference", "sequences", "logs"]:
        (out / d).mkdir(exist_ok=True)
    status = dict(
        status="running",
        package_version=__version__,
        methods=list(methods),
        subset=subset,
        discovery_parameters="native defaults",
        context_bp=2000,
        flank_bp=500,
    )
    write_json(out / "run.json", status)
    runner = Runner(out, bins)
    try:
        discover(samples, out, runner, methods)
        features, maps, nproteins = prepare_annotations(samples, out)
        families = collect(samples, out, methods, features, maps)
        if not discovery_only:
            annotate(families, out, runner, nproteins, threads)
            genomes = {s["id"]: fasta_dict(out / s["fasta"]) for s in samples}
            totals = finish_boundaries(out, runner, families, genomes, features, subset)
            status["summary"] = totals
        status["status"] = "complete"
        status["stage"] = "native_only" if discovery_only else "boundaries"
        # Native outputs are immutable throughout all post-processing.
        for rel, expected in read_json(out / "native_manifest.json")[-1][
            "sha256"
        ].items():
            if digest(out / rel) != expected:
                raise ValueError("Native output was modified: " + rel)
    except Exception as error:
        status.update(status="failed", error=str(error))
        raise
    finally:
        write_json(out / "run.json", status)
    return status
