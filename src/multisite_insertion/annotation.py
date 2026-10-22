"""Post-hoc deposited annotations and translated protein evidence; no discovery gates."""

import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .io import completed_blast, fasta_dict, write_json

BLASTX_FORMAT = "7 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qframe stitle"


def prepare_annotations(samples, out):
    out = Path(out)
    features, mappings, proteins = {}, {}, []
    protein_ids = set()
    for sample in samples:
        sid = sample["id"]
        genome = fasta_dict(out / sample["fasta"])
        aliases = defaultdict(set)
        for contig in genome:
            for alias in {contig, *contig.split("|")}:
                aliases[alias].add(contig)
        features[sid], mappings[sid] = [], {}
        for file in sample["annotations"]:
            with (out / file).open() as handle:
                annotated_records = list(
                    SeqIO.parse(handle, sample.get("annotation_format", "embl"))
                )
            for record in annotated_records:
                candidates = aliases[record.id]
                if len(candidates) != 1:
                    raise ValueError(
                        f"Annotation ID {record.id} must match one FASTA ID or pipe-separated accession in {sid}"
                    )
                contig = next(iter(candidates))
                if contig in mappings[sid]:
                    raise ValueError(f"Duplicate annotation for {sid}/{contig}")
                if str(record.seq).upper() != str(genome[contig].seq).upper():
                    raise ValueError(
                        f"Annotation/FASTA sequence mismatch: {sid}/{contig}"
                    )
                mappings[sid][contig] = record.id
                for i, feature in enumerate(record.features):
                    if feature.type not in [
                        "CDS",
                        "rRNA",
                        "tRNA",
                        "mobile_element",
                        "repeat_region",
                    ]:
                        continue
                    if feature.location is None:
                        raise ValueError(
                            f"Annotation has no usable location: {record.id}, feature {i}"
                        )
                    rid = f"{record.id}_feature_{i}"
                    row = dict(
                        id=rid,
                        accession=record.id,
                        type=feature.type,
                        start0=int(feature.location.start),
                        end0=int(feature.location.end),
                        parts=[
                            dict(start0=int(part.start), end0=int(part.end))
                            for part in feature.location.parts
                        ],
                        feature_length_bp=len(feature.location),
                        strand=feature.location.strand,
                        location=str(feature.location),
                        qualifiers={
                            k: v
                            for k, v in feature.qualifiers.items()
                            if k != "translation"
                        },
                    )
                    features[sid].append(row)
                    if "translation" in feature.qualifiers:
                        pid = sid + "__" + rid
                        if pid in protein_ids:
                            raise ValueError("Duplicate protein ID: " + pid)
                        protein_ids.add(pid)
                        proteins.append(
                            SeqRecord(
                                Seq(feature.qualifiers["translation"][0]),
                                id=pid,
                                description=" ".join(
                                    feature.qualifiers.get("product", ["unknown"])
                                ),
                            )
                        )
        write_json(out / "reference" / f"{sid}_features.json", features[sid])
    (out / "reference").mkdir(exist_ok=True)
    SeqIO.write(proteins, out / "reference/deposited_proteins.faa", "fasta")
    write_json(out / "reference/accession_maps.json", mappings)
    return features, mappings, len(proteins)


def annotate(families, out, runner, protein_count, threads=2):
    out = Path(out).resolve()
    hits_file = out / "results/deposited.blastx.tsv"
    query = out / "inputs/all_native_families.fasta"
    query_count = sum(r["consensus_available"] for r in families)
    if protein_count and query_count:
        db = out / "reference/deposited"
        runner.run(
            "deposited_makeblastdb",
            "makeblastdb",
            [
                "-in",
                out / "reference/deposited_proteins.faa",
                "-dbtype",
                "prot",
                "-parse_seqids",
                "-out",
                db,
            ],
        )
        runner.run(
            "deposited_blastx",
            "blastx",
            [
                "-query",
                query,
                "-db",
                db,
                "-evalue",
                "1e-5",
                "-max_target_seqs",
                "1000",
                "-num_threads",
                threads,
                "-outfmt",
                BLASTX_FORMAT,
                "-out",
                hits_file,
            ],
        )
        completed_blast(hits_file, query_count)
    else:
        hits_file.write_text(
            "# Not run: no translated annotation proteins or no family consensuses.\n"
        )
    hits = defaultdict(list)
    for line in hits_file.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        hits[f[0]].append(
            dict(
                subject=f[1],
                identity_pct=float(f[2]),
                alignment_aa=int(f[3]),
                query_start=int(f[4]),
                query_end=int(f[5]),
                evalue=float(f[8]),
                bitscore=float(f[9]),
                subject_length=int(f[11]),
                title=f[13],
            )
        )
    for row in families:
        hs = sorted(hits[row["query"]], key=lambda h: -h["bitscore"])
        tnp = [h for h in hs if re.search("transposase", h["title"], re.I)]
        row["best_deposited_protein_hit"] = hs[0] if hs else None
        row["best_deposited_tnp_hit"] = tnp[0] if tnp else None
        row["native_records_overlapping_tnp_CDS"] = sum(
            any(
                a["type"] == "CDS"
                and re.search("transposase", " ".join(a["product"]), re.I)
                for a in o["annotations"]
            )
            for o in row["occurrences"]
        )
        row["tnp_supported"] = bool(tnp and row["native_records_overlapping_tnp_CDS"])
        row["deposited_tnp_products"] = sorted(
            {
                " ".join(a["product"])
                for o in row["occurrences"]
                for a in o["annotations"]
                if a["type"] == "CDS"
                and re.search("transposase", " ".join(a["product"]), re.I)
            }
        )
    write_json(out / "results/all_families_annotated.json", families)
    return families
