"""S2: catalog, joint FNA index, measured distances and the S3 pair table for one verified pool.

catalog  One row per BioSample from pool.tsv + S1 verify.json + the BioSample attributes in the S0 NCBI snapshot:
         FNA path/SHA256, records, plasmids, taxon, linked read types, annotation status, GFF feature coordinates
         outside their contig, and the batch50 source rule with two fixes (SOURCE_RULES; source_class is a flag,
         nothing is removed from discovery). Writes the joint FASTA (headers sample|assembly|contig, each contig kept
         separate, bases unchanged) and contig_index.tsv.gz (to --index-dir when given: large pools keep it on scratch,
         joint_index.json records its path and SHA256). Pool rows whose primary download is not verified are left out
         of the catalog and listed in catalog_excluded.tsv with their S1 status and reason.
sources  Recompute the source flags in catalog.tsv only (no joint FNA rewrite).
ani      skani triangle -E --slow over every FNA (same command as batch50); edge list kept on scratch.
pairs    Clone groups (ANI >= 99.99 and AF >= 95 both ways), species clusters (single linkage, ANI >= 95 and
         AF >= 50 both ways), genus from the organism name. Pair table for S3: per assembly up to 1 clone partner,
         up to 3 nearest non-clone same-species neighbours, and the nearest other-species partner with skani output;
         plus audit sets: all pairs of one single-project and one multi-project cluster (10-20 members), and
         same-genus cross-project pairs in ANI tiers including pairs with no skani output. Each pair records whether
         it meets the species threshold itself (the cluster label is single linkage). Samples in no pair are
         no_pair_selected_pending_review, not "cannot be compared". These selection numbers are
         initial settings for the first pool, not calibrated thresholds.
"""

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "callout_scale_v1" / "batch50"))
import build_panel as bp  # noqa: E402  (batch50 source rules, with the fixes below)
from build_panel import SKANI, attrs, components  # noqa: E402

# Source rules v2 = batch50 rules with two fixes for this pool: plant keywords need word boundaries
# ("stem" matched "system"), and facility phrases such as "wastewater treatment plant" are not plants.
# Patched only in this process; batch50's frozen panel is unchanged (none of its plant labels used these words).
SOURCE_RULES = "batch50_source_class+v2(word-bounded plant keywords; '<facility> plant' -> facility)"
bp.PLANT_TEXT = re.compile(r"\broot|rhizosph|nodule|\bleaf\b|leaves|xylem|\bplants?\b|\bseeds?\b|\bstems?\b|potato|tomato|"
                           r"banana|\brice\b|wheat|maize|citri|vetiver|juniperus|\barum\b|\bfruits?\b|ginger|turmeric", re.I)
FACILITY = re.compile(r"\b(treatment|sewage|refinery|power|processing|wastewater|waste ?water|industrial|chemical|"
                      r"manufacturing|desalination|water)[- ]plant\b", re.I)

CLONE_ANI, CLONE_AF = 99.99, 95.0
SPECIES_ANI, SPECIES_AF = 95.0, 50.0
N_SPECIES_NEIGHBOURS = 3
AUDIT_MIN, AUDIT_MAX = 10, 20
NO_RUN = dict(host="", host_scientific_name="", isolation_source="", environmental_medium="")
DISTANCE_TIERS = [("ani_99.99_100", 99.99, 100.01), ("ani_99_99.99", 99.0, 99.99), ("ani_97_99", 97.0, 99.0),
                  ("ani_95_97", 95.0, 97.0), ("ani_90_95", 90.0, 95.0), ("ani_lt90", 0.0, 90.0)]
PER_TIER = 5


def source_class(report):
    """Source flag from BioSample attributes (raw fields are kept in the catalog)."""
    at = {k: FACILITY.sub(lambda m: m.group(1) + " facility", v) for k, v in attrs(report).items()}
    return bp.source_class(at, NO_RUN, report["assembly_info"]["biosample"].get("package", ""))


def load_reports(work):
    reports = {}
    for src in ("refseq", "genbank"):
        for line in open(str(work / "raw" / ("ncbi_%s.jsonl" % src))):
            d = json.loads(line)
            reports[d["accession"]] = d
    return reports


def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_tsv(path, rows):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "wt") as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def gff_out_of_range(path, lengths):
    """GFF3 feature rows whose seqid is not an FNA record or whose coordinates fall outside it."""
    n = 0
    with open(str(path)) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split("\t", 5)
            L = lengths.get(c[0])
            if L is None or not 1 <= int(c[3]) <= int(c[4]) <= L:
                n += 1
    return n


def catalog(a):
    data, s2 = a.work / "data", a.work / "s2"
    s2.mkdir(exist_ok=True)
    pool = read_tsv(a.pool_dir / "pool.tsv")
    reports = load_reports(a.work)
    state = {}
    for line in open(str(data / "state.jsonl")):
        d = json.loads(line)
        state[d["accession"]] = d
    rows, contigs, excluded = [], [], []
    joint = s2 / "joint.fna.gz"
    with gzip.open(str(joint), "wt", compresslevel=6) as out:
        for p in pool:
            acc = p["download_accession"]
            vp = data / "genomes" / acc / "verify.json"
            v = json.loads(vp.read_text()) if vp.exists() else {}
            if state.get(acc, {}).get("status") != "verified" or v.get("status") != "verified":
                excluded.append(dict(sample_id=p["biosample"], assembly=acc, s1_status=state.get(acc, {}).get("status", "pending"),
                                     reason=state.get(acc, {}).get("reason", "") or v.get("reason", ""),
                                     note="primary download not verified; not used by any caller"))
                continue
            fna = [f["file"] for f in v["files"] if f["file"].endswith("_genomic.fna")][0]
            rep = reports.get(acc) or reports.get(p["genbank_accession"])
            at = attrs(rep)
            cls, note = source_class(rep)
            sample = p["biosample"]
            first = len(contigs)
            with open(str(data / "genomes" / acc / fna)) as f:
                name, seq = None, []
                for line in itertools.chain(f, [">"]):
                    if line.startswith(">"):
                        if name is not None:
                            s = "".join(seq)
                            header = "%s|%s|%s" % (sample, acc, name)
                            out.write(">" + header + "\n")
                            for i in range(0, len(s), 80):
                                out.write(s[i:i + 80] + "\n")
                            contigs.append(dict(sample_id=sample, assembly=acc, contig=name, length=len(s),
                                                sha256=hashlib.sha256(s.upper().encode()).hexdigest(), joint_header=header))
                        name, seq = (line[1:].split() or [None])[0], []
                    else:
                        seq.append(line.strip())
            mol = v.get("molecule_types", {})
            lengths = {c["contig"]: c["length"] for c in contigs[first:]}
            gff = data / "genomes" / acc / "genomic.gff"
            oor = gff_out_of_range(gff, lengths) if gff.exists() else ""
            rows.append(dict(
                sample_id=sample, assembly=acc, genbank_accession=p["genbank_accession"], role=p["role"],
                organism=p["organism"], genus=p["organism"].split()[0].strip("[]"), tax_id=p["tax_id"], strain=p["strain"],
                bioproject=p["bioproject"], assembly_level=p["assembly_level"], records=v["records"], total_bp=v["fna_bp"],
                plasmid_records=mol.get("Plasmid", 0), annotation=v.get("annotation", ""),
                annotation_status=v.get("annotation_status", ""), gff_features_out_of_range=oor, linked_read_types=p["linked_read_types"],
                n_usable_runs=p["n_usable_runs"], atypical=p["atypical"], source_class=cls, source_note=note, source_rules=SOURCE_RULES,
                host=at.get("host", ""), isolation_source=at.get("isolation_source", ""), env_medium=at.get("env_medium", ""),
                geo_loc_name=at.get("geo_loc_name", ""), collection_date=at.get("collection_date", ""),
                checkm_completeness=rep.get("checkm_info", {}).get("completeness", ""),
                checkm_contamination=rep.get("checkm_info", {}).get("contamination", ""),
                fna_path="genomes/%s/%s" % (acc, fna), fna_sha256=[f["sha256"] for f in v["files"] if f["file"] == fna][0]))
    # Each joint record must equal its source contig.
    by = {(c["assembly"], c["contig"]): c for c in contigs}
    for r in rows:
        v = json.loads((data / "genomes" / r["assembly"] / "verify.json").read_text())
        for k, h in v["record_sha256"].items():
            assert by[(r["assembly"], k)]["sha256"] == h, (r["assembly"], k)
    write_tsv(a.out / "catalog.tsv", rows)
    index = (a.index_dir or a.out) / "contig_index.tsv.gz"
    index.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(index, contigs)
    if excluded:
        write_tsv(a.out / "catalog_excluded.tsv", excluded)
    (s2 / "genomes.list").write_text("".join(str(data / r["fna_path"]) + "\n" for r in rows))
    meta = dict(joint_fna=str(joint), joint_fna_sha256=sha256(joint), samples=len(rows), contigs=len(contigs),
                contig_bp=sum(c["length"] for c in contigs), header="sample_id|assembly|contig",
                check="every joint record SHA256 equals the S1 record SHA256 of its source FNA",
                pool_rows=len(pool), excluded_not_verified=len(excluded),
                annotation_status=dict(Counter(r["annotation_status"] for r in rows)),
                gff_features_out_of_range=sum(r["gff_features_out_of_range"] or 0 for r in rows))
    if a.index_dir:
        meta.update(contig_index=str(index), contig_index_sha256=sha256(index))
    (a.out / "joint_index.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2), dict(Counter(r["source_class"] for r in rows)))


def sources(a):
    """Recompute source flags in catalog.tsv from the S0 snapshot without rewriting the joint FNA."""
    reports = load_reports(a.work)
    cat = read_tsv(a.out / "catalog.tsv")
    changed = Counter()
    for r in cat:
        rep = reports.get(r["assembly"]) or reports.get(r["genbank_accession"])
        cls, note = source_class(rep)
        if cls != r["source_class"]:
            changed["%s->%s" % (r["source_class"], cls)] += 1
        r.update(source_class=cls, source_note=note, source_rules=SOURCE_RULES)
    cols = list(cat[0])
    for c in ("source_rules",):
        if c in cols:
            cols.remove(c)
    cols.insert(cols.index("source_note") + 1, "source_rules")
    write_tsv(a.out / "catalog.tsv", [{c: r[c] for c in cols} for r in cat])
    print(dict(Counter(r["source_class"] for r in cat)), "changed:", dict(changed))


def ani(a, threads):
    s2 = a.work / "s2"
    edges = s2 / "skani_edges.tsv"
    cmd = [SKANI, "triangle", "-l", str(s2 / "genomes.list"), "-E", "--slow", "-t", str(threads), "-o", str(edges)]
    subprocess.run(cmd, check=True)
    ver = subprocess.run([SKANI, "--version"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    meta = dict(command=" ".join(cmd), version=ver, binary_sha256=sha256(SKANI), edges=str(edges),
                edges_sha256=sha256(edges), rows=sum(1 for _ in open(str(edges))) - 1)
    (a.out / "skani_run.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


def is_linked(r):
    """Has existing ONT/PacBio run links (the baseline role "linked" is exactly this; assembly-first pools have one role)."""
    return bool(r["linked_read_types"])


def stable(key):
    return hashlib.sha256(key.encode()).hexdigest()


def pairs(a):
    cat = read_tsv(a.out / "catalog.tsv")
    path_to = {str(a.work / "data" / r["fna_path"]): r["sample_id"] for r in cat}
    ids = [r["sample_id"] for r in cat]
    info = {r["sample_id"]: r for r in cat}
    E = {}
    for e in csv.DictReader(open(str(a.work / "s2" / "skani_edges.tsv")), delimiter="\t"):
        x, y = path_to[e["Ref_file"]], path_to[e["Query_file"]]
        E[frozenset((x, y))] = dict(ani=float(e["ANI"]), af={x: float(e["Align_fraction_ref"]), y: float(e["Align_fraction_query"])})

    def meets(e, ani, af):
        return e is not None and e["ani"] >= ani and min(e["af"].values()) >= af

    clone = components(ids, [tuple(k) for k, v in E.items() if meets(v, CLONE_ANI, CLONE_AF)])
    # Single linkage: members of one cluster need not all meet the threshold pairwise (see pair_meets_species_threshold).
    species = components(ids, [tuple(k) for k, v in E.items() if meets(v, SPECIES_ANI, SPECIES_AF)])

    chosen = defaultdict(set)

    def add(x, y, why):
        chosen[frozenset((x, y))].add(why)

    nb = defaultdict(list)
    for k, v in E.items():
        x, y = tuple(k)
        nb[x].append((v["ani"], min(v["af"].values()), y))
        nb[y].append((v["ani"], min(v["af"].values()), x))
    for x in ids:
        ranked = sorted(nb[x], reverse=True)
        cl = [y for _, _, y in ranked if clone[y] == clone[x]]
        if cl:
            add(x, cl[0], "clone_partner")
        for y in [y for _, _, y in ranked if species[y] == species[x] and clone[y] != clone[x]][:N_SPECIES_NEIGHBOURS]:
            add(x, y, "species_neighbour")
        other = [y for _, _, y in ranked if species[y] != species[x]]
        if other:
            add(x, other[0], "nearest_other_species")

    # Audit sets. (1) one-project cluster: program consistency; (2) cross-project cluster: neighbour-selection misses;
    # (3) same-genus distance tiers from different projects, including pairs with no skani output, as controls.
    sizes = Counter(species.values())
    members = defaultdict(list)
    for i in ids:
        members[species[i]].append(i)
    nproj = {c: len({info[i]["bioproject"] for i in m}) for c, m in members.items()}
    linked = {c: sum(is_linked(info[i]) for i in m) for c, m in members.items()}
    eligible = [c for c, n in sizes.items() if AUDIT_MIN <= n <= AUDIT_MAX]
    single = sorted((c for c in eligible if nproj[c] == 1), key=lambda c: (-linked[c], -sizes[c], c))
    multi = sorted((c for c in eligible if nproj[c] > 1), key=lambda c: (-nproj[c], -linked[c], c))
    audit = dict(single_project=single[0] if single else "", multi_project=multi[0] if multi else "")
    for tag, c in audit.items():
        if c:
            for x, y in itertools.combinations(sorted(members[c]), 2):
                add(x, y, "audit_%s_all_vs_all" % tag)
    by_genus = defaultdict(list)
    for i in ids:
        by_genus[info[i]["genus"]].append(i)
    tier_pool = defaultdict(list)
    for g, m in by_genus.items():
        for x, y in itertools.combinations(sorted(m), 2):
            if info[x]["bioproject"] == info[y]["bioproject"]:
                continue
            e = E.get(frozenset((x, y)))
            if e is None:
                tier = "no_skani_output_same_genus"
            else:
                tier = [t for t, lo, hi in DISTANCE_TIERS if lo <= e["ani"] < hi][0]
            tier_pool[tier].append((x, y))
    for tier, ps in tier_pool.items():
        ps.sort(key=lambda p: (not (is_linked(info[p[0]]) or is_linked(info[p[1]])), stable("|".join(p))))
        for x, y in ps[:PER_TIER]:
            add(x, y, "distance_tier_" + tier)

    rows = []
    for k, why in sorted(chosen.items(), key=lambda kv: sorted(kv[0])):
        x, y = sorted(k)
        e = E.get(k)
        rows.append(dict(
            pair_id="P%05d" % (len(rows) + 1), sample_a=x, sample_b=y, assembly_a=info[x]["assembly"], assembly_b=info[y]["assembly"],
            reasons=";".join(sorted(why)), ani=e["ani"] if e else "", af_a=e["af"][x] if e else "", af_b=e["af"][y] if e else "",
            skani="hit" if e else "no_output", same_clone=clone[x] == clone[y],
            same_species_cluster=species[x] == species[y],
            pair_meets_species_threshold=meets(e, SPECIES_ANI, SPECIES_AF),
            same_genus_name=info[x]["genus"] == info[y]["genus"], same_bioproject=info[x]["bioproject"] == info[y]["bioproject"],
            both_linked=is_linked(info[x]) and is_linked(info[y]),
            either_linked=is_linked(info[x]) or is_linked(info[y])))
    in_pair = Counter(s for r in rows for s in (r["sample_a"], r["sample_b"]))
    for r in cat:
        sid = r["sample_id"]
        r.update(clone_group=clone[sid], species_cluster=species[sid], pairs_selected=in_pair[sid],
                 pair_status="in_pair_table" if in_pair[sid] else "no_pair_selected_pending_review",
                 skani_output_pairs=len(nb[sid]))
    write_tsv(a.out / "catalog.tsv", cat)
    write_tsv(a.out / "pairs.tsv", rows)
    write_tsv(a.out / "groups.tsv", [dict(sample_id=r["sample_id"], clone_group=r["clone_group"], species_cluster=r["species_cluster"],
                                          genus=r["genus"], role=r["role"], bioproject=r["bioproject"]) for r in cat])
    same_cluster_below = [r["pair_id"] for r in rows if r["same_species_cluster"] and not r["pair_meets_species_threshold"]]
    summ = dict(
        schema="fna_pool_s2_v2", samples=len(ids), skani_pairs_reported=len(E), possible_pairs=len(ids) * (len(ids) - 1) // 2,
        samples_without_skani_output=[i for i in ids if not nb[i]],
        samples_without_skani_output_note=("no skani output under its default screen (~80% ANI k-mer screen) and AF >= 15 reporting; "
                                           "not evidence that no local homology exists"),
        pair_status=dict(Counter(r["pair_status"] for r in cat)),
        pair_status_linked=dict(Counter(r["pair_status"] for r in cat if is_linked(r))),
        clone_groups=len(set(clone.values())), clone_groups_size_gt1=sum(1 for n in Counter(clone.values()).values() if n > 1),
        species_clusters=len(sizes), species_cluster_sizes=dict(sorted(Counter(sizes.values()).items())),
        species_clusters_ge_10=sum(1 for n in sizes.values() if n >= 10), genera=len(by_genus),
        source_rules=SOURCE_RULES, source_class=dict(Counter(r["source_class"] for r in cat)),
        source_class_linked=dict(Counter(r["source_class"] for r in cat if is_linked(r))),
        pairs=len(rows), pairs_with_multiple_reasons=sum(";" in r["reasons"] for r in rows),
        pairs_by_reason=dict(sorted(Counter(x for r in rows for x in r["reasons"].split(";")).items())),
        pairs_same_cluster_below_species_threshold=same_cluster_below,
        pairs_both_linked=sum(r["both_linked"] for r in rows), pairs_either_linked=sum(r["either_linked"] for r in rows),
        audit_clusters={tag: dict(id=c, size=sizes.get(c, 0), linked=linked.get(c, 0), bioprojects=nproj.get(c, 0),
                                  organisms=dict(Counter(info[i]["organism"] for i in members.get(c, []))))
                        for tag, c in audit.items()},
        distance_tier_candidates={t: len(v) for t, v in sorted(tier_pool.items())},
        settings=dict(clone="ANI >= %.2f and AF >= %.0f both" % (CLONE_ANI, CLONE_AF),
                      species="single linkage over pairs with ANI >= %.0f and AF >= %.0f both; a cluster label does not mean "
                              "every member pair meets it (pair_meets_species_threshold)" % (SPECIES_ANI, SPECIES_AF),
                      neighbours="1 clone partner, %d nearest same-cluster non-clone, 1 nearest other-cluster (skani output)"
                                 % N_SPECIES_NEIGHBOURS,
                      audit="all pairs of one single-project and one multi-project cluster with %d-%d members" % (AUDIT_MIN, AUDIT_MAX),
                      distance_tiers="same genus, different BioProject; up to %d pairs per tier %s + no skani output; "
                                     "linked first, then SHA256 order" % (PER_TIER, [t for t, _, _ in DISTANCE_TIERS]),
                      note="initial settings for this pool, not calibrated thresholds; S3 must judge each pair from its own alignment"))
    (a.out / "s2_summary.json").write_text(json.dumps(summ, indent=2) + "\n")
    print(json.dumps(summ, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["catalog", "sources", "ani", "pairs"])
    p.add_argument("--work", type=Path, required=True, help="pool work dir on scratch (raw/, data/, s2/)")
    p.add_argument("--pool-dir", type=Path, required=True, help="repo dir with pool.tsv")
    p.add_argument("--out", type=Path, required=True, help="repo dir for S2 tables")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--index-dir", type=Path, help="catalog: write contig_index.tsv.gz here (scratch) instead of --out")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    dict(catalog=catalog, sources=sources, pairs=pairs)[a.mode](a) if a.mode != "ani" else ani(a, a.threads)


if __name__ == "__main__":
    main()
