"""Per-sample and per-family call-out report for the five-sample acceptance run (task step 6).

All native records are kept. Layers (length, boundary, flanks, unambiguous triplet,
v1 and audited read support, tnp status) are columns, never filters. Exact
coordinate aliases and partial overlaps are counted separately; support from
different tools at one place is never added.
"""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

PRIORITY_BOUNDARIES = ("deposited_mobile_element", "provisional_homology_median")


def fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def tnp_status(fam):
    if not fam.get("consensus_available", True):
        return "unavailable"
    return "supported" if fam.get("tnp_supported") else "unsupported"


def report(sample, callout, v1_evidence, audited):
    res = callout / "results"
    families = {f["query"]: f for f in json.loads((res / "all_families_annotated.json").read_text())}
    native_records = sum(f["raw_occurrences"] for f in families.values())
    members = json.loads((res / "member_boundaries.json").read_text())
    invalid = json.loads((res / "invalid_native_records.json").read_text())
    seqs = {k: fasta(callout / "sequences" / (k + ".fasta"))
            for k in ("candidate_insertion", "upstream_500", "downstream_500", "triplet")}
    v1 = {r["member"]: r for r in csv.DictReader(open(v1_evidence / "member_read_support.tsv"), delimiter="\t")}
    aud = {r["member"]: r for r in csv.DictReader(open(audited / "member_audit.tsv"), delimiter="\t")}

    coords = Counter((m["accession"], m["candidate_start0"], m["candidate_end0"]) for m in members)
    by_contig = defaultdict(list)
    for m in members:
        by_contig[m["accession"]].append(m)
    rows = []
    for m in members:
        key = (m["accession"], m["candidate_start0"], m["candidate_end0"])
        partial = sum(1 for o in by_contig[m["accession"]]
                      if o is not m and (o["accession"], o["candidate_start0"], o["candidate_end0"]) != key
                      and o["candidate_start0"] < m["candidate_end0"] and o["candidate_end0"] > m["candidate_start0"])
        length = m["candidate_end0"] - m["candidate_start0"]
        trip = seqs["triplet"].get(m["member"], "")
        missing = [k for k, v in seqs.items() if m["member"] not in v]
        a = aud.get(m["member"], {})
        v = v1.get(m["member"], {})
        rows.append(dict(
            sample=sample, member=m["member"], family=m["family"], method=m["method"], contig=m["accession"],
            start0=m["candidate_start0"], end0=m["candidate_end0"], length=length, strand=m["strand"],
            length_layer=">=400" if length >= 400 else "<400", boundary_source=m["boundary_source"],
            flanks_complete=bool(m["flanks_complete"]), triplet_acgt_only=bool(trip) and set(trip.upper()) <= set("ACGT"),
            exported_fasta_missing=",".join(missing), exact_coordinate_aliases=coords[key] - 1,
            partial_overlap_records=partial, tnp=tnp_status(families[m["family"]]),
            v1_class=v.get("class", "missing"), v1_full=v.get("full", ""),
            v1_full_ge2=v.get("full", "0").isdigit() and int(v.get("full", 0)) >= 2,
            audited_class=a.get("audited_class", "missing"), audited_full_units=a.get("units_full_support", ""),
            audited_full=a.get("audited_class", "").startswith("audited_full_span"),
            priority=length >= 400 and bool(m["flanks_complete"]) and m["boundary_source"] in PRIORITY_BOUNDARIES
            and bool(trip) and set(trip.upper()) <= set("ACGT"),
        ))
    fam_rows = []
    for fam, rs in sorted(defaultdict(list, {f: [r for r in rows if r["family"] == f] for f in families}).items()):
        f = families[fam]
        fam_rows.append(dict(
            sample=sample, family=fam, method=f["method"], tnp=tnp_status(f), native_records=f["raw_occurrences"],
            exported_records=len(rs), long_records=sum(r["length_layer"] == ">=400" for r in rs),
            exact_coordinate_groups=len({(r["contig"], r["start0"], r["end0"]) for r in rs}),
            flanks_complete=sum(r["flanks_complete"] for r in rs), triplet_acgt_only=sum(r["triplet_acgt_only"] for r in rs),
            priority_records=sum(r["priority"] for r in rs), v1_full_ge2=sum(r["v1_full_ge2"] for r in rs),
            audited_full=sum(r["audited_full"] for r in rs),
            boundary_sources=";".join("%s=%d" % kv for kv in sorted(Counter(r["boundary_source"] for r in rs).items())),
        ))
    layer = lambda pred: dict(records=sum(1 for r in rows if pred(r)),
                              v1_full_ge2=sum(1 for r in rows if pred(r) and r["v1_full_ge2"]),
                              audited_full=sum(1 for r in rows if pred(r) and r["audited_full"]))
    summary = dict(
        sample=sample, native_families=len(families), native_records=native_records, exported_records=len(rows),
        invalid_native_records=len(invalid),
        exact_coordinate_groups=len(coords), records_with_exact_aliases=sum(r["exact_coordinate_aliases"] > 0 for r in rows),
        records_with_partial_overlap=sum(r["partial_overlap_records"] > 0 for r in rows),
        boundary_sources=dict(Counter(r["boundary_source"] for r in rows)),
        flanks_complete=sum(r["flanks_complete"] for r in rows), triplet_acgt_only=sum(r["triplet_acgt_only"] for r in rows),
        exported_fasta_missing=sum(bool(r["exported_fasta_missing"]) for r in rows),
        tnp_records=dict(Counter(r["tnp"] for r in rows)),
        layers={"all": layer(lambda r: True), "<400": layer(lambda r: r["length_layer"] == "<400"),
                ">=400": layer(lambda r: r["length_layer"] == ">=400"), "priority": layer(lambda r: r["priority"])},
        v1_classes=dict(Counter(r["v1_class"] for r in rows)), audited_classes=dict(Counter(r["audited_class"] for r in rows)),
        v1_full_ge2_but_not_audited=sum(r["v1_full_ge2"] and not r["audited_full"] for r in rows),
        audited_but_not_v1=sum(r["audited_full"] and not r["v1_full_ge2"] for r in rows),
        note="Records, not insertion events. Exact aliases and partial overlaps are diagnostics, not locus counts.",
    )
    return rows, fam_rows, summary


def write_tsv(path, rows):
    with open(path, "w") as h:
        if not rows:
            return
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, default=Path(__file__).resolve().parent / "execution_manifest.json")
    p.add_argument("--proj", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout"))
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    all_rows, all_fams, summaries = [], [], []
    for e in json.load(open(a.manifest))["samples"]:
        rows, fams, summary = report(e["sample_id"], Path(e["v1_callout"]), Path(e["v1_evidence"]),
                                     a.proj / "callout_scale_v1" / e["sample_id"] / "evidence_audited")
        all_rows += rows
        all_fams += fams
        summaries.append(summary)
    write_tsv(a.out / "members.tsv", all_rows)
    write_tsv(a.out / "families.tsv", all_fams)
    (a.out / "samples.json").write_text(json.dumps(summaries, indent=2) + "\n")
    for s in summaries:
        print(s["sample"], s["exported_records"], s["layers"]["priority"], s["audited_classes"])


if __name__ == "__main__":
    main()
