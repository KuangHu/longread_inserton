"""Coordinate, orientation, retention and immutable native-output checks."""

from pathlib import Path

from .io import digest, fasta_dict, read_json, samples_in_run, write_json


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(out):
    out = Path(out).resolve()
    require(read_json(out / "run.json")["status"] == "complete", "Run did not complete")
    samples = samples_in_run(out)
    for path, expected in read_json(out / "native_manifest.json")[-1]["sha256"].items():
        require(digest(out / path) == expected, "Native output changed: " + path)
    rows = read_json(out / "results/member_boundaries.json")
    families = read_json(out / "results/selected_families.json")
    invalid = read_json(out / "results/invalid_native_records.json")
    require(
        len(rows) + len(invalid) == sum(f["raw_occurrences"] for f in families),
        "Selected record count mismatch",
    )
    expected_ids = {
        f'{f["query"]}_m{i}'
        for f in families
        for i in range(1, f["raw_occurrences"] + 1)
    }
    require(
        {r["member"] for r in rows + invalid} == expected_ids,
        "Selected member IDs changed",
    )
    genome = {s["id"]: fasta_dict(out / s["fasta"]) for s in samples}
    originals = {
        f'{f["query"]}_m{i}': (f["organism"], o)
        for f in families
        for i, o in enumerate(f["occurrences"], 1)
    }
    kinds = [
        "candidate_insertion",
        "upstream_500",
        "downstream_500",
        "triplet",
        "genomic_core",
        "genomic_left_500",
        "genomic_right_500",
    ]
    seqs = {k: fasta_dict(out / "sequences" / (k + ".fasta")) for k in kinds}
    for records in seqs.values():
        require(
            set(records) == {r["member"] for r in rows}, "Export member IDs changed"
        )
    for r in rows:
        q = r["member"]
        sid, o = originals[q]
        require(
            (r["native_start0"], r["native_end0"], r["strand"])
            == (o["start0"], o["end0"], o["strand"]),
            q + ": native interval changed",
        )
        g = genome[sid][o["contig"]].seq
        a, b = r["candidate_start0"], r["candidate_end0"]
        require(0 <= a < b <= len(g), q + ": candidate out of bounds")
        core = g[a:b]
        require(
            str(seqs["genomic_core"][q].seq) == str(core), q + ": genomic core mismatch"
        )
        require(
            str(seqs["genomic_left_500"][q].seq) == str(g[max(0, a - 500) : a]),
            q + ": left mismatch",
        )
        require(
            str(seqs["genomic_right_500"][q].seq) == str(g[b : b + 500]),
            q + ": right mismatch",
        )
        if r["strand"] == "-":
            core = core.reverse_complement()
        require(
            str(seqs["candidate_insertion"][q].seq) == str(core),
            q + ": oriented core mismatch",
        )
        for tag, field in [
            ("upstream_500", "upstream"),
            ("downstream_500", "downstream"),
        ]:
            flank = g[r[field + "_genomic_start0"] : r[field + "_genomic_end0"]]
            if r["strand"] == "-":
                flank = flank.reverse_complement()
            require(str(seqs[tag][q].seq) == str(flank), q + ": flank mismatch")
        require(
            str(seqs["triplet"][q].seq)
            == str(seqs["upstream_500"][q].seq + core + seqs["downstream_500"][q].seq),
            q + ": triplet mismatch",
        )
    result = dict(
        status="pass",
        selected_families=len(families),
        selected_native_records=len(rows) + len(invalid),
        exported_records=len(rows),
        invalid_records_retained=len(invalid),
        native_outputs_unchanged=True,
        scope="Data consistency, not biological accuracy or transposition activity",
    )
    write_json(out / "results/validation.json", result)
    return result
