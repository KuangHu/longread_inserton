"""Pilot boundary estimator, retained without a new biological candidate filter.

Discovery intervals, inferred homology extents and reference annotations remain
separate. No empty-site, TIR or TSD requirement is introduced.
"""

from pathlib import Path
from collections import defaultdict
from statistics import median_low
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json, re
from .io import fasta_dict


def core_link_count(qa, sa, btop, qmeta, smeta):
    qp, sp = qa, sa
    n = 0
    for token in re.findall(r"\d+|..", btop):
        if token.isdigit():
            k = int(token)
            lo = max(0, qmeta["core_start0"] - qp, smeta["core_start0"] - sp)
            hi = min(k, qmeta["core_end0"] - qp, smeta["core_end0"] - sp)
            n += max(0, hi - lo)
            qp += k
            sp += k
        else:
            qc, sc = token
            if (
                qc != "-"
                and sc != "-"
                and qmeta["core_start0"] <= qp < qmeta["core_end0"]
                and smeta["core_start0"] <= sp < smeta["core_end0"]
            ):
                n += 1
            qp += qc != "-"
            sp += sc != "-"
    return n, qp, sp


def infer(selected, genomes, features, p, meta):
    p = Path(p)
    contexts = fasta_dict(p / "inputs/selected_member_contexts.fasta")
    best = {}
    raw_counts = defaultdict(int)
    for l in (p / "results/member_contexts_btop.blastn.tsv").read_text().splitlines():
        if l.startswith("#"):
            continue
        f = l.split("\t")
        q, s = f[:2]
        if (
            q not in meta
            or s not in meta
            or q == s
            or meta[q]["family"] != meta[s]["family"]
        ):
            continue
        raw_counts[q] += 1
        qa, qb = int(f[4]) - 1, int(f[5])
        sa, sb = int(f[6]) - 1, int(f[7])
        if qa >= qb or sa >= sb:
            continue
        linked, qp, sp = core_link_count(qa, sa, f[12], meta[q], meta[s])
        assert qp == qb and sp == sb, (q, s, qp, qb, sp, sb)
        if linked == 0:
            continue
        row = dict(
            query=q,
            subject=s,
            start0=qa,
            end0=qb,
            subject_start0=sa,
            subject_end0=sb,
            identity_pct=float(f[2]),
            alignment_columns=int(f[3]),
            evalue=float(f[8]),
            bitscore=float(f[9]),
            core_to_core_aligned_bases=linked,
            btop=f[12],
        )
        if (q, s) not in best or (linked, row["bitscore"]) > (
            best[q, s]["core_to_core_aligned_bases"],
            best[q, s]["bitscore"],
        ):
            best[q, s] = row
    observations = defaultdict(list)
    for r in best.values():
        observations[r["query"]].append(r)
    records = []
    seqout = {
        name: []
        for name in [
            "upstream_500",
            "candidate_insertion",
            "downstream_500",
            "triplet",
            "flanks_for_blast",
            "genomic_left_500",
            "genomic_right_500",
            "genomic_core",
        ]
    }
    for q, m in meta.items():
        fam = selected[m["family"]]
        idx = int(q.rsplit("_m", 1)[1]) - 1
        o = fam["occurrences"][idx]
        assert o["contig"] == m["contig"] and o["start0"] == m["native_start0"]
        hs = observations[q]
        seq = genomes[m["organism"]][m["contig"]].seq
        strand = o["strand"]
        row = dict(
            member=q,
            family=m["family"],
            method=fam["method"],
            organism=m["organism"],
            accession=o["accession"],
            strand=strand,
            native_start0=o["start0"],
            native_end0=o["end0"],
            native_length=o["end0"] - o["start0"],
            peer_count=len(hs),
            available_peers=fam["raw_occurrences"] - 1,
            homology_HSPs=hs,
        )
        if hs:
            starts = [h["start0"] for h in hs]
            ends = [h["end0"] for h in hs]
            lo = median_low(starts)
            hi = median_low(ends)
            a, b = (
                (m["extract_start0"] + lo, m["extract_start0"] + hi)
                if strand == "+"
                else (m["extract_end0"] - hi, m["extract_end0"] - lo)
            )
            row.update(
                inferred_start0=a,
                inferred_end0=b,
                inferred_length=b - a,
                oriented_start_spread_bp=max(starts) - min(starts),
                oriented_end_spread_bp=max(ends) - min(ends),
                oriented_start_observed_min0=min(starts),
                oriented_start_observed_max0=max(starts),
                oriented_end_observed_min0=min(ends),
                oriented_end_observed_max0=max(ends),
                context_window_touched=any(
                    h["start0"] == 0 or h["end0"] == len(contexts[q]) for h in hs
                ),
                median_boundary_touches_window=lo == 0 or hi == len(contexts[q]),
            )
        else:
            a, b = o["start0"], o["end0"]
            row.update(inferred_start0=None, inferred_end0=None, inferred_length=None)
        # Retain all overlapping deposited full-element annotations for review.
        overlaps = []
        for feat in features[m["organism"]]:
            if feat["type"] != "mobile_element" or feat["accession"] != o["accession"]:
                continue
            ov = sum(
                max(0, min(o["end0"], part["end0"]) - max(o["start0"], part["start0"]))
                for part in feat["parts"]
            )
            if ov > 0:
                overlaps.append(
                    dict(
                        id=feat["id"],
                        start0=feat["start0"],
                        end0=feat["end0"],
                        name=feat["qualifiers"].get("mobile_element_type", []),
                        overlap_native_bp=ov,
                        contiguous=len(feat["parts"]) == 1,
                    )
                )
        row["overlapping_mobile_annotations"] = overlaps
        # Anchor an annotation to the native interval midpoint, rather than adopting
        # a neighboring element that shares only a few terminal bases. Keep every overlap.
        midpoint = (o["start0"] + o["end0"] - 1) // 2
        anchored = [
            ann
            for ann in overlaps
            if ann["contiguous"] and ann["start0"] <= midpoint < ann["end0"]
        ]
        row["midpoint_anchored_mobile_annotations"] = anchored
        if len(anchored) == 1:
            ann = anchored[0]
            fa, fb = ann["start0"], ann["end0"]
            source = "deposited_mobile_element"
            if hs:
                row.update(
                    inferred_minus_annotation_start_bp=a - fa,
                    inferred_minus_annotation_end_bp=b - fb,
                )
        elif hs:
            fa, fb = a, b
            source = (
                "homology_extent_censored_at_window"
                if row["median_boundary_touches_window"]
                else "provisional_homology_median"
            )
        else:
            fa, fb = a, b
            source = "native_interval_only_unresolved"
        row.update(
            candidate_start0=fa,
            candidate_end0=fb,
            candidate_start1=fa + 1,
            candidate_end1=fb,
            candidate_length=fb - fa,
            boundary_source=source,
        )
        assert 0 <= fa < fb <= len(seq)
        la = max(0, fa - 500)
        rb = min(len(seq), fb + 500)
        left = seq[la:fa]
        core = seq[fa:fb]
        right = seq[fb:rb]
        if strand == "+":
            up, ins, down = left, core, right
            upcoords = (la, fa)
            downcoords = (fb, rb)
        else:
            up, ins, down = (
                right.reverse_complement(),
                core.reverse_complement(),
                left.reverse_complement(),
            )
            upcoords = (fb, rb)
            downcoords = (la, fa)
        row.update(
            upstream_genomic_start0=upcoords[0],
            upstream_genomic_end0=upcoords[1],
            downstream_genomic_start0=downcoords[0],
            downstream_genomic_end0=downcoords[1],
            upstream_length=len(up),
            downstream_length=len(down),
            flanks_complete=len(up) == 500 and len(down) == 500,
        )
        # Describe exact short repeats only; no TSD/TIR requirement or boundary search.
        row["exact_flank_direct_repeat_lengths"] = [
            k
            for k in range(2, 21)
            if len(up) >= k and len(down) >= k and str(up[-k:]) == str(down[:k])
        ]
        row["exact_terminal_inverted_repeat_lengths"] = [
            k
            for k in range(4, min(51, len(ins) // 2 + 1))
            if str(ins[:k]) == str(ins[-k:].reverse_complement())
        ]
        row["candidate_overlapping_tnp_CDS"] = []
        for feat in features[m["organism"]]:
            if (
                feat["type"] != "CDS"
                or feat["accession"] != o["accession"]
                or not re.search(
                    "transposase", " ".join(feat["qualifiers"].get("product", [])), re.I
                )
            ):
                continue
            ov = sum(
                max(0, min(fb, part["end0"]) - max(fa, part["start0"]))
                for part in feat["parts"]
            )
            if ov > 0:
                row["candidate_overlapping_tnp_CDS"].append(
                    dict(
                        id=feat["id"],
                        product=feat["qualifiers"].get("product", []),
                        overlap_bp=ov,
                        CDS_length_bp=feat["feature_length_bp"],
                        fully_contained=ov == feat["feature_length_bp"],
                    )
                )
        row["review_notes"] = []
        if overlaps and not anchored:
            row["review_notes"].append(
                "No unique contiguous mobile-element annotation contains the native interval midpoint; overlapping annotation boundaries were not adopted."
            )
        if not hs:
            row["review_notes"].append(
                "No peer HSP links both native cores; no homology boundary was inferred."
            )
        if source == "native_interval_only_unresolved":
            row["review_notes"].append(
                "The exported core is the original native interval only."
            )
        if not row["flanks_complete"]:
            row["review_notes"].append(
                "Flanks are clipped at a contig end; circular origins are not wrapped."
            )
        if row.get("context_window_touched"):
            row["review_notes"].append(
                "At least one peer HSP reaches the analysis window edge; its outer extent is censored."
            )
        if source == "homology_extent_censored_at_window":
            row["review_notes"].append(
                "The exported homology extent touches the window edge and is not a complete insertion boundary."
            )
        desc = f"candidate={o['accession']}:{fa+1}-{fb} family_orientation={strand} boundary={source}"
        for kind, s in [
            ("upstream_500", up),
            ("candidate_insertion", ins),
            ("downstream_500", down),
            ("triplet", up + ins + down),
        ]:
            seqout[kind].append(SeqRecord(s, id=q, description=desc))
        for kind, s, coords in [
            ("genomic_left_500", left, (la, fa)),
            ("genomic_right_500", right, (fb, rb)),
            ("genomic_core", core, (fa, fb)),
        ]:
            seqout[kind].append(
                SeqRecord(
                    s,
                    id=q,
                    description=f"segment={o['accession']}:{coords[0]+1}-{coords[1]} genome_forward=+ boundary={source}",
                )
            )
        seqout["flanks_for_blast"].extend(
            [
                SeqRecord(up, id=q + "__up", description=""),
                SeqRecord(down, id=q + "__down", description=""),
            ]
        )
        records.append(row)
    for name, rr in seqout.items():
        SeqIO.write(rr, p / "sequences" / f"{name}.fasta", "fasta")
    SeqIO.write(
        [contexts[q] for q in meta],
        p / "inputs/selected_member_contexts.fasta",
        "fasta",
    )
    (p / "results/member_boundaries.json").write_text(
        json.dumps(records, indent=2) + "\n"
    )
    (p / "results/selected_families.json").write_text(
        json.dumps(list(selected.values()), indent=2) + "\n"
    )
    return records
