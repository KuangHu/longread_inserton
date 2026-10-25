"""S3 small-scale pilot: choose pairs across distances, projects and selection reasons; pick loci to plot; summarise.

select   From s2/pairs.tsv: one pair per distance tier (same genus, different BioProject), one clone partner, one
         cross-project species neighbour with both samples linked, one nearest-other-species pair, one pair from each
         audit cluster. Linked pairs first, then pair_id. Writes s3_pilot/pilot_pairs.tsv with FNA/GFF3 paths.
plotlist Loci to draw: per pair up to 2 anchored cores < 2500 bp, 1 anchored >= 2500 bp, 1 replacement, 1 unresolved.
summary  Per-pair and overall counts from calls/loci/one-flank tables.
allpairs Every pair of s2/pairs.tsv in the pilot_pairs.tsv format (pilot_reason = "all"), for the full S3 run.
crosstab Short anchored loci (< SHORT bp): CDS content of the core from the carrier GFF3 (1-2 complete CDS with no
         CDS/rRNA/tRNA crossing either boundary, 0 complete CDS, > 2, or boundary_cuts_feature), joined with
         difference_class, boundary consistency and (optional) read-check evidence level. Writes short_loci_review.tsv.
         This is a review table, not S4: candidate cores still need boundary/repeat review.
"""

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

PICK = [("distance_tier_ani_99.99_100", None), ("distance_tier_ani_99_99.99", None), ("distance_tier_ani_97_99", None),
        ("distance_tier_ani_95_97", None), ("distance_tier_ani_90_95", None), ("distance_tier_ani_lt90", None),
        ("distance_tier_no_skani_output_same_genus", None), ("clone_partner", "linked"),
        ("species_neighbour", "linked_cross_project"), ("nearest_other_species", "linked"),
        ("audit_multi_project_all_vs_all", None), ("audit_single_project_all_vs_all", None)]
SHORT = 2500


def read(path):
    return list(csv.DictReader(open(str(path)), delimiter="\t"))


def write(path, rows):
    with open(str(path), "w") as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def pair_row(p, cat, data, reason):
    row = dict(pair_id=p["pair_id"], pilot_reason=reason, reasons=p["reasons"], ani=p["ani"], af_a=p["af_a"], af_b=p["af_b"],
               same_bioproject=p["same_bioproject"], both_linked=p["both_linked"])
    for k, s in (("a", p["sample_a"]), ("b", p["sample_b"])):
        c = cat[s]
        d = Path(c["fna_path"]).parent
        row.update({"sample_%s" % k: s, "organism_%s" % k: c["organism"], "bioproject_%s" % k: c["bioproject"],
                    "level_%s" % k: c["assembly_level"], "fna_%s" % k: str(data / c["fna_path"]),
                    "gff_%s" % k: str(data / d / "genomic.gff")})
    return row


def allpairs(a):
    pairs, cat = read(a.pool_dir / "s2" / "pairs.tsv"), {r["sample_id"]: r for r in read(a.pool_dir / "s2" / "catalog.tsv")}
    out = [pair_row(p, cat, a.data, "all") for p in pairs]
    a.out.mkdir(parents=True, exist_ok=True)
    write(a.out / "pairs.tsv", out)
    print(len(out), "pairs;", len({r["sample_a"] for r in out} | {r["sample_b"] for r in out}), "samples")


def crosstab(a):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from s3_pair_calls import read_gff_cds
    gff = {}
    for p in read(a.pilot):
        gff[p["sample_a"]], gff[p["sample_b"]] = p["gff_a"], p["gff_b"]
    ev = {l["locus_id"]: l for l in read(a.readcheck)} if a.readcheck else {}
    feats, rows = {}, []
    for l in read(a.loci):
        if l["status"] != "two_flank_anchored" or int(l["core_len"]) >= SHORT:
            continue
        if l["carrier"] not in feats:
            feats[l["carrier"]] = read_gff_cds(gff[l["carrier"]])
        s, e = int(l["core_start"]), int(l["core_end"])
        fs = feats[l["carrier"]][l["contig"]]
        full = [f for f in fs if f[3] == "CDS" and s <= f[0] and f[1] <= e]
        cut = [f for f in fs if f[0] < s < f[1] or f[0] < e < f[1]]
        cds = ("boundary_cuts_feature" if cut else "0_complete_CDS" if not full else
               "1-2_complete_CDS" if len(full) <= 2 else ">2_complete_CDS")
        r = ev.get(l["locus_id"], {})
        rows.append(dict(locus_id=l["locus_id"], carrier=l["carrier"], contig=l["contig"], core_start=s, core_end=e,
                         core_len=l["core_len"], difference_class=l["difference_class"], boundary_consistency=l["boundary_consistency"],
                         anchored_support=l["anchored_support"], tandem_like=l["tandem_like"], cds_class=cds,
                         complete_cds=len(full), pseudo_cds=sum(1 for f in full if f[5]),
                         cut_features=";".join("%s:%s" % (f[3], f[4]) for f in cut), products=";".join(f[4] for f in full),
                         evidence_level=r.get("evidence_level", ""), carrier_count_unit=r.get("carrier_count_unit", "")))
    write(a.out / "short_loci_review.tsv", rows)
    tab = Counter((r["evidence_level"], r["difference_class"], r["cds_class"]) for r in rows)
    for k, v in sorted(tab.items()):
        print("\t".join(k), v, sep="\t")


def select(a):
    pairs, cat = read(a.pool_dir / "s2" / "pairs.tsv"), {r["sample_id"]: r for r in read(a.pool_dir / "s2" / "catalog.tsv")}
    out, used = [], set()
    for reason, extra in PICK:
        cand = [p for p in pairs if reason in p["reasons"].split(";") and p["pair_id"] not in used]
        if extra == "linked":
            cand = [p for p in cand if p["both_linked"] == "True"] or cand
        if extra == "linked_cross_project":
            cand = [p for p in cand if p["both_linked"] == "True" and p["same_bioproject"] == "False"] or cand
        cand.sort(key=lambda p: (p["either_linked"] != "True", p["pair_id"]))
        if not cand:
            continue
        p = cand[0]
        used.add(p["pair_id"])
        out.append(pair_row(p, cat, a.data, reason))
    a.out.mkdir(parents=True, exist_ok=True)
    write(a.out / "pilot_pairs.tsv", out)
    for r in out:
        print(r["pair_id"], r["pilot_reason"], r["ani"] or "no_skani", r["organism_a"], "|", r["organism_b"], r["same_bioproject"])


def plotlist(a):
    loci, calls = read(a.loci), read(a.calls)
    first_pair = {}
    for c in calls:
        first_pair.setdefault(c["locus_id"], c["pair_id"])
    by = defaultdict(list)
    for l in loci:
        by[first_pair[l["locus_id"]]].append(l)
    pick = []
    for pid, ls in sorted(by.items()):
        ls.sort(key=lambda l: l["locus_id"])
        anchored = [l for l in ls if l["status"] == "two_flank_anchored"]
        pick += [l for l in anchored if int(l["core_len"]) < SHORT][:2]
        pick += [l for l in anchored if int(l["core_len"]) >= SHORT][:1]
        pick += [l for l in ls if l["status"] == "replacement"][:1]
        pick += [l for l in ls if l["status"] == "unresolved"][:1]
    print(" ".join(l["locus_id"] for l in pick))


def summary(a):
    pilot = read(a.pilot)
    loci, calls = read(a.loci), read(a.calls)
    one = []
    for p in pilot:
        f = a.work / p["pair_id"] / "one_flank.tsv"
        if f.exists():
            one += read(f) if f.stat().st_size > 20 else []
    lid = {l["locus_id"]: l for l in loci}
    pair_loci, pair_calls = defaultdict(set), Counter()
    for c in calls:
        pair_loci[c["pair_id"]].add(c["locus_id"])
        pair_calls[c["pair_id"]] += 1
    one_by = Counter(o["pair_id"] for o in one)
    rows = []
    for p in pilot:
        pid = p["pair_id"]
        ls = [lid[x] for x in sorted(pair_loci[pid])]
        an = [l for l in ls if l["status"] == "two_flank_anchored"]
        rows.append(dict(pair_id=pid, pilot_reason=p["pilot_reason"], reasons=p["reasons"], ani=p["ani"] or "no_skani",
                         same_bioproject=p["same_bioproject"], organism_a=p["organism_a"], organism_b=p["organism_b"],
                         calls=pair_calls[pid], loci=len(ls), anchored=len(an),
                         anchored_lt2500=sum(1 for l in an if int(l["core_len"]) < SHORT),
                         anchored_tandem=sum(1 for l in an if l["tandem_like"] == "True"),
                         replacement=sum(1 for l in ls if l["status"] == "replacement"),
                         unresolved=sum(1 for l in ls if l["status"] == "unresolved"),
                         one_flank=one_by[pid]))
    write(a.out / "pilot_summary.tsv", rows)
    reasons = Counter(r for l in loci if l["status"] == "unresolved" for r in l["reasons"].split(";") if r)
    tot = dict(pairs=len(pilot), calls=len(calls), loci=len(loci), status=dict(Counter(l["status"] for l in loci)),
               anchored_lt2500=sum(1 for l in loci if l["status"] == "two_flank_anchored" and int(l["core_len"]) < SHORT),
               loci_seen_in_multiple_pairs=sum(1 for l in loci if int(l["n_pairs"]) > 1),
               unresolved_reasons=dict(reasons.most_common()), one_flank=len(one),
               core_len_anchored=sorted(int(l["core_len"]) for l in loci if l["status"] == "two_flank_anchored"))
    (a.out / "pilot_summary.json").write_text(json.dumps(tot, indent=1) + "\n")
    for r in rows:
        print("\t".join(str(v) for v in r.values()))
    print(json.dumps({k: v for k, v in tot.items() if k != "core_len_anchored"}, indent=1))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["select", "plotlist", "summary", "allpairs", "crosstab"])
    p.add_argument("--readcheck", type=Path, help="readcheck_loci.tsv (crosstab)")
    p.add_argument("--pool-dir", type=Path)
    p.add_argument("--data", type=Path)
    p.add_argument("--work", type=Path)
    p.add_argument("--pilot", type=Path)
    p.add_argument("--loci", type=Path)
    p.add_argument("--calls", type=Path)
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    dict(select=select, plotlist=plotlist, summary=summary, allpairs=allpairs, crosstab=crosstab)[a.mode](a)


if __name__ == "__main__":
    main()
