"""Descriptive flank and cross-tool summaries; no candidate rejection."""

from pathlib import Path
from collections import defaultdict, Counter
from itertools import combinations
import json


def summarize(p):
    p = Path(p)
    rows = json.loads((p / "results/member_boundaries.json").read_text())
    byid = {r["member"]: r for r in rows}
    families = defaultdict(list)
    for r in rows:
        families[r["family"]].append(r)
    hits = defaultdict(list)
    for line in (p / "results/flanks_all_vs_all.blastn.tsv").read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        q, qs = f[0].rsplit("__", 1)
        s, ss = f[1].rsplit("__", 1)
        if q == s or qs != ss or byid[q]["family"] != byid[s]["family"]:
            continue
        # Reciprocal BLAST observations are retained as such in the raw output.
        if q > s:
            continue
        h = dict(
            identity_pct=float(f[2]),
            alignment_columns=int(f[3]),
            query_span_bp=abs(int(f[5]) - int(f[4])) + 1,
            subject_span_bp=abs(int(f[7]) - int(f[6])) + 1,
            evalue=float(f[8]),
            bitscore=float(f[9]),
            query_start1=int(f[4]),
            query_end1=int(f[5]),
            subject_start1=int(f[6]),
            subject_end1=int(f[7]),
        )
        hits[q, s, qs].append(h)
    pairs = []
    summaries = []
    for fam, rr in families.items():
        fp = []
        for aa, bb in combinations(rr, 2):
            a, b = sorted([aa, bb], key=lambda r: r["member"])
            pair = dict(family=fam, member_a=a["member"], member_b=b["member"])
            for side, longside in [("up", "upstream"), ("down", "downstream")]:
                hh = hits[a["member"], b["member"], side]
                best = (
                    max(
                        hh,
                        key=lambda h: (
                            min(h["query_span_bp"], h["subject_span_bp"]),
                            h["bitscore"],
                        ),
                    )
                    if hh
                    else None
                )
                ov = (
                    max(
                        0,
                        min(
                            a[longside + "_genomic_end0"], b[longside + "_genomic_end0"]
                        )
                        - max(
                            a[longside + "_genomic_start0"],
                            b[longside + "_genomic_start0"],
                        ),
                    )
                    if a["organism"] == b["organism"]
                    and a["accession"] == b["accession"]
                    else 0
                )
                pair[side] = dict(
                    HSP_count=len(hh),
                    longest_shared_span_HSP=best,
                    physical_window_overlap_bp=ov,
                    has_HSP_spanning_at_least_250bp_in_both=any(
                        min(h["query_span_bp"], h["subject_span_bp"]) >= 250 for h in hh
                    ),
                )
            fp.append(pair)
            pairs.append(pair)
        ss = dict(
            family=fam,
            organism=rr[0]["organism"],
            method=rr[0]["method"],
            native_records=len(rr),
            candidate_length_min=min(r["candidate_length"] for r in rr),
            candidate_length_max=max(r["candidate_length"] for r in rr),
            boundary_sources=dict(Counter(r["boundary_source"] for r in rr)),
            records_containing_full_annotated_tnp=sum(
                any(t["fully_contained"] for t in r["candidate_overlapping_tnp_CDS"])
                for r in rr
            ),
            record_pairs=len(fp),
            pairs_no_corresponding_flank_HSP=sum(
                x["up"]["HSP_count"] == 0 and x["down"]["HSP_count"] == 0 for x in fp
            ),
            pairs_neither_flank_HSP_spans_250bp=sum(
                not x["up"]["has_HSP_spanning_at_least_250bp_in_both"]
                and not x["down"]["has_HSP_spanning_at_least_250bp_in_both"]
                for x in fp
            ),
            pairs_physical_flank_overlap=sum(
                x["up"]["physical_window_overlap_bp"] > 0
                or x["down"]["physical_window_overlap_bp"] > 0
                for x in fp
            ),
            max_oriented_start_spread_bp=max(
                r.get("oriented_start_spread_bp", 0) for r in rr
            ),
            max_oriented_end_spread_bp=max(
                r.get("oriented_end_spread_bp", 0) for r in rr
            ),
        )
        summaries.append(ss)
    cross = []
    for r in rows:
        if r["method"] != "RS":
            continue
        for q in rows:
            if (
                q["method"] != "PILER"
                or q["organism"] != r["organism"]
                or q["accession"] != r["accession"]
            ):
                continue
            ov = max(
                0,
                min(q["native_end0"], r["native_end0"])
                - max(q["native_start0"], r["native_start0"]),
            )
            if ov:
                cross.append(
                    dict(
                        RS_member=r["member"],
                        PILER_member=q["member"],
                        native_overlap_bp=ov,
                        RS_minus_PILER_candidate_start_bp=r["candidate_start0"]
                        - q["candidate_start0"],
                        RS_minus_PILER_candidate_end_bp=r["candidate_end0"]
                        - q["candidate_end0"],
                        same_family_orientation=r["strand"] == q["strand"],
                    )
                )
    for name, obj in [
        ("flank_pair_comparison", pairs),
        ("family_summary", summaries),
        ("cross_tool_locus_overlaps", cross),
    ]:
        (p / "results" / f"{name}.json").write_text(json.dumps(obj, indent=2) + "\n")
    cols = [
        "member",
        "family",
        "accession",
        "strand",
        "native_start0",
        "native_end0",
        "inferred_start0",
        "inferred_end0",
        "candidate_start1",
        "candidate_end1",
        "candidate_length",
        "boundary_source",
        "peer_count",
        "upstream_genomic_start0",
        "upstream_genomic_end0",
        "downstream_genomic_start0",
        "downstream_genomic_end0",
    ]
    (p / "results/member_coordinates.txt").write_text(
        "\t".join(cols)
        + "\n"
        + "\n".join("\t".join(str(r.get(k, "")) for k in cols) for r in rows)
        + "\n"
    )
    (p / "results/candidates.bed").write_text(
        "".join(
            f"{r['accession']}\t{r['candidate_start0']}\t{r['candidate_end0']}\t{r['member']}\t0\t{r['strand']}\n"
            for r in rows
        )
    )
    cols = [
        "family",
        "native_records",
        "candidate_length_min",
        "candidate_length_max",
        "records_containing_full_annotated_tnp",
        "record_pairs",
        "pairs_no_corresponding_flank_HSP",
        "pairs_neither_flank_HSP_spans_250bp",
        "pairs_physical_flank_overlap",
        "max_oriented_start_spread_bp",
        "max_oriented_end_spread_bp",
    ]
    (p / "results/family_summary.txt").write_text(
        "\t".join(cols)
        + "\n"
        + "\n".join("\t".join(str(r[k]) for k in cols) for r in summaries)
        + "\n"
    )
    return summaries
