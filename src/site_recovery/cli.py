"""Explicit-input, fresh-output long-read insertion recovery workflow."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pysam

from .core import Config, candidates_from_vcf, recover, write_outputs


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_reference(bam, ref):
    expected = dict(zip(ref.references, ref.lengths))
    observed = dict(zip(bam.references, bam.lengths))
    if expected != observed:
        raise ValueError("BAM reference names/lengths must exactly match the supplied assembly")
    for sq in bam.header.to_dict().get("SQ", []):
        if "M5" in sq:
            md5 = hashlib.md5(ref.fetch(sq["SN"]).upper().encode()).hexdigest()
            if sq["M5"].lower() != md5:
                raise ValueError(f"BAM reference MD5 mismatch: {sq['SN']}")
    if not bam.has_index():
        raise ValueError("BAM must be indexed and coordinate sorted")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["run", "recover"])
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--platform", choices=["hifi", "ccs", "ont"], required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--reads", type=Path, nargs="+")
    p.add_argument("--bam", type=Path)
    p.add_argument("--vcf", type=Path)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--sniffles", default="sniffles")
    p.add_argument("--mosaic-af-min", type=float, default=0.05)
    for name in ("flank", "anchor", "min_mapq", "min_support", "min_other_support",
                 "min_size", "max_size", "search_slop", "max_sequences", "max_reads", "repeat_limit"):
        p.add_argument("--"+name.replace("_", "-"), type=int, default=getattr(Config(), name))
    p.add_argument("--length-tolerance", type=float, default=Config().length_tolerance)
    return p


def execute(args):
    names = set(Config.__dataclass_fields__) - {"platform"}
    cfg = Config(platform="hifi" if args.platform == "ccs" else args.platform,
                 **{k: getattr(args, k) for k in names})
    cfg.validate()
    if args.threads < 1 or not 0 < args.mosaic_af_min < 1:
        raise ValueError("threads must be positive and mosaic-af-min must be in (0,1)")
    if args.mode == "run":
        if not args.reads or args.bam or args.vcf:
            raise ValueError("run requires --reads and does not accept --bam/--vcf")
        input_paths = [args.reference] + args.reads
    else:
        if not args.bam or not args.vcf or args.reads:
            raise ValueError("recover requires --bam and --vcf and does not accept --reads")
        input_paths = [args.reference, args.bam, args.vcf]
    for path in input_paths:
        if not path.is_file():
            raise ValueError(f"Input does not exist: {path}")
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be new or empty; stale output reuse is disabled")
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"status": "running", "sample": args.sample, "mode": args.mode,
                "config": asdict(cfg), "commands": [], "inputs": [],
                "software": {k: importlib.metadata.version(k)
                             for k in ("insertion-site-recovery", "pysam", "edlib", "pyspoa")}}

    def save():
        (out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")

    def command(cmd, log, stdout=None):
        manifest["commands"].append(cmd)
        save()
        with (out/log).open("w") as err:
            subprocess.run(cmd, check=True, stdout=stdout if stdout is not None else err, stderr=err)

    save()
    try:
        for path in input_paths:
            manifest["inputs"].append({"path": str(path.resolve()), "sha256": sha256_file(path)})
        reference_hash = manifest["inputs"][0]["sha256"]
        # Stage a private reference so indexing never changes the supplied inputs.
        reference = out/"reference.fasta"
        if args.reference.suffix.lower() in {".gz", ".bgz"}:
            raise ValueError("Supply uncompressed reference FASTA; compressed reads are supported")
        shutil.copyfile(args.reference, reference)
        pysam.faidx(str(reference))
        if args.mode == "run":
            for tool in (args.minimap2, args.sniffles):
                if not shutil.which(tool):
                    raise ValueError(f"Executable not found: {tool}")
            for key, tool in (("minimap2", args.minimap2), ("sniffles", args.sniffles)):
                manifest["software"][key] = subprocess.check_output(
                    [tool, "--version"], text=True, stderr=subprocess.STDOUT).strip()
            preset = "map-hifi" if cfg.platform == "hifi" else "map-ont"
            sam = out/"alignments.unsorted.sam"
            bam_path, vcf_path = out/"alignments.bam", out/"calls.vcf"
            with sam.open("w") as fh:
                command([args.minimap2, "-a", "-x", preset, "-t", str(args.threads),
                         str(reference)] + [str(p.resolve()) for p in args.reads], "minimap2.log", fh)
            pysam.sort("-@", str(args.threads), "-o", str(bam_path), str(sam))
            pysam.index(str(bam_path))
            sam.unlink()
            command([args.sniffles, "--input", str(bam_path), "--vcf", str(vcf_path),
                     "--reference", str(reference), "--threads", str(args.threads),
                     "--all-contigs", "--mosaic", "--mosaic-include-germline",
                     "--mosaic-af-min", str(args.mosaic_af_min), "--max-svlen-mosaic", str(cfg.max_size),
                     "--minsupport", str(cfg.min_support), "--minsvlen", str(cfg.min_size),
                     "--mapq", str(cfg.min_mapq), "--output-rnames"], "sniffles.log")
            manifest["reference_validation"] = "alignment generated here against staged, hashed reference"
        else:
            bam_path, vcf_path = args.bam.resolve(), args.vcf.resolve()
            manifest["reference_validation"] = "header names/lengths and available M5; absent M5 cannot prove sequence identity"
        with pysam.FastaFile(str(reference)) as ref, pysam.AlignmentFile(str(bam_path), "rb") as bam:
            validate_reference(bam, ref)
            records = [recover(c, bam, ref, cfg, reference_hash, args.sample)
                       for c in candidates_from_vcf(vcf_path)]
        summary = write_outputs(records, out)
        manifest.update(status="complete", summary=summary,
                        warning="Synthetic-tested prototype; local allele support does not establish strain or mobility")
        save()
        return summary
    except Exception as exc:
        manifest.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        save()
        raise


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        print(json.dumps(execute(args), ensure_ascii=False))
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
