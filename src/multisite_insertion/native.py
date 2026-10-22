"""Direct native discovery and lossless coordinate normalization."""

import re
import shutil
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from .io import digest, fasta_dict, write_json


def discover(samples, out, runner, methods):
    out = Path(out).resolve()
    manifest = []
    for s in samples:
        for method in methods:
            dest = out / "native" / s["id"] / method
            source = s.get("native", {}).get(method)
            if source:
                shutil.copytree(source, dest)
                manifest.append(
                    {
                        "sample": s["id"],
                        "method": method,
                        "source": source,
                        "mode": "imported",
                    }
                )
            else:
                dest.mkdir(parents=True)
                seq = out / s["fasta"]
                prefix = s["id"] + "_" + method
                if method == "RS":
                    runner.run(
                        prefix + "_lmer",
                        "build_lmer_table",
                        ["-sequence", seq, "-freq", "lmers.freq"],
                        dest,
                    )
                    runner.run(
                        prefix + "_discover",
                        "RepeatScout",
                        [
                            "-sequence",
                            seq,
                            "-freq",
                            "lmers.freq",
                            "-output",
                            "families.fasta",
                            "-ranges",
                            "ranges.tsv",
                        ],
                        dest,
                    )
                else:
                    runner.run(
                        prefix + "_pals",
                        "pals",
                        ["-self", seq, "-out", "hits.gff", "-log", "pals.log"],
                        dest,
                    )
                    runner.run(
                        prefix + "_piler",
                        "piler2",
                        [
                            "-trs",
                            "hits.gff",
                            "-out",
                            "families.gff",
                            "-piles",
                            "piles.gff",
                            "-images",
                            "images.gff",
                            "-log",
                            "piler.log",
                        ],
                        dest,
                    )
                    for name in ["fams", "aligned_fams", "cons"]:
                        (dest / name).mkdir()
                    # Empty native discovery is valid, not an execution failure.
                    if any(
                        l and not l.startswith("#")
                        for l in (dest / "families.gff").read_text().splitlines()
                    ):
                        runner.run(
                            prefix + "_extract",
                            "piler2",
                            [
                                "-trs2fasta",
                                "families.gff",
                                "-seq",
                                seq,
                                "-path",
                                "fams",
                            ],
                            dest,
                        )
                        for f in sorted((dest / "fams").iterdir()):
                            if not f.is_file():
                                continue
                            label = f.name
                            runner.run(
                                prefix + "_muscle_" + label,
                                "muscle",
                                [
                                    "-in",
                                    f,
                                    "-out",
                                    dest / "aligned_fams" / label,
                                    "-maxiters",
                                    "1",
                                    "-diags1",
                                ],
                                dest,
                            )
                            runner.run(
                                prefix + "_cons_" + label,
                                "piler2",
                                [
                                    "-cons",
                                    dest / "aligned_fams" / label,
                                    "-out",
                                    dest / "cons" / label,
                                    "-label",
                                    label,
                                ],
                                dest,
                            )
                    with (dest / "library.fasta").open("wb") as handle:
                        for f in sorted((dest / "cons").iterdir()):
                            if f.is_file():
                                handle.write(f.read_bytes())
                manifest.append(
                    {"sample": s["id"], "method": method, "mode": "native_defaults"}
                )
            required = (
                ["library.fasta", "families.gff"]
                if method == "PILER"
                else ["families.fasta", "ranges.tsv"]
            )
            for name in required:
                if not (dest / name).is_file():
                    raise ValueError(f"Missing native output {dest / name}")
    manifest.append(
        {
            "sha256": {
                str(p.relative_to(out)): digest(p)
                for p in sorted((out / "native").rglob("*"))
                if p.is_file()
            }
        }
    )
    write_json(out / "native_manifest.json", manifest)


def collect(samples, out, methods, features, accession_maps):
    families, sequences = [], []
    for s in samples:
        genome = fasta_dict(Path(out) / s["fasta"])
        for method in methods:
            folder = Path(out) / "native" / s["id"] / method
            occurrences = defaultdict(list)
            file = folder / ("families.gff" if method == "PILER" else "ranges.tsv")
            for number, line in enumerate(file.read_text().splitlines(), 1):
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if method == "PILER":
                    match = re.search(r"Family ([0-9.]+)", fields[8])
                    if not match:
                        raise ValueError(
                            f"Unknown PILER family syntax in {file}:{number}"
                        )
                    family = match.group(1)
                    a, b, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
                else:
                    family = fields[3]
                    a, b, strand = int(fields[1]), int(fields[2]), fields[4]
                contig = fields[0]
                if contig not in genome:
                    raise ValueError(
                        f'Native contig {contig} is absent from {s["id"]} FASTA'
                    )
                if strand not in ("+", "-"):
                    raise ValueError(f"Unexpected native strand: {strand}")
                accession = accession_maps[s["id"]].get(contig, contig)
                annotations = []
                for feature in features[s["id"]]:
                    if feature["accession"] != accession:
                        continue
                    overlap = overlap_parts(a, b, feature["parts"])
                    if overlap:
                        annotations.append(
                            dict(
                                feature_id=feature["id"],
                                type=feature["type"],
                                overlap_bp=overlap,
                                feature_length_bp=feature["feature_length_bp"],
                                product=feature["qualifiers"].get("product", []),
                                locus_tag=feature["qualifiers"].get("locus_tag", []),
                            )
                        )
                occurrences[family].append(
                    dict(
                        contig=contig,
                        accession=accession,
                        start0=a,
                        end0=b,
                        strand=strand,
                        in_bounds=0 <= a < b <= len(genome[contig]),
                        annotations=annotations,
                    )
                )
            consensus_path = folder / (
                "library.fasta" if method == "PILER" else "families.fasta"
            )
            consensus = fasta_dict(consensus_path)
            # Preserve families with no consensus as well as consensus-only families.
            for label in dict.fromkeys([*consensus, *occurrences]):
                safe = label.replace("=", "")
                if not re.fullmatch(r"[A-Za-z0-9_.-]+", safe):
                    raise ValueError("Unexpected native family label: " + label)
                q = f'{s["id"]}_{method}_{safe}'
                record = consensus.get(label)
                if record is not None:
                    record.id = q
                    record.description = ""
                    sequences.append(record)
                oo = occurrences[label]
                families.append(
                    dict(
                        query=q,
                        organism=s["id"],
                        method=method,
                        native_family=label,
                        length=len(record) if record is not None else None,
                        raw_occurrences=len(oo),
                        occurrences=oo,
                        consensus_available=record is not None,
                    )
                )
    if len({r["query"] for r in families}) != len(families):
        raise ValueError("Normalized family IDs collide")
    SeqIO.write(sequences, Path(out) / "inputs/all_native_families.fasta", "fasta")
    write_json(Path(out) / "results/all_native_families.json", families)
    return families


def overlap_parts(a, b, parts):
    """Compound-origin genes must never use their chromosome-wide bounding box."""
    return sum(max(0, min(b, part["end0"]) - max(a, part["start0"])) for part in parts)
