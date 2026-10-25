"""S2: catalog, joint FNA index, measured distances and the S3 pair table for one verified pool.

catalog  One row per BioSample from pool.tsv + S1 verify.json + the BioSample attributes in the S0 NCBI snapshot:
         FNA path/SHA256, records, plasmids, taxon, linked read types, and the batch50 source rule
         (source_class is a flag here; nothing is removed from discovery). Writes the joint FASTA
         (headers sample|assembly|contig, each contig kept separate, bases unchanged) and contig_index.tsv.
ani      skani triangle -E --slow over every FNA (same command as batch50); edge list kept on scratch.
pairs    Clone groups (ANI >= 99.99 and AF >= 95 both ways), species clusters (single linkage, ANI >= 95 and
         AF >= 50 both ways), genus from the organism name. Pair table for S3: per assembly up to 1 clone partner,
         up to 3 nearest non-clone same-species neighbours, and the nearest other-species partner with skani output;
         plus all pairs of one audit cluster (10-20 members) to measure what neighbour selection misses.
         Pairs without skani output are recorded as no_hit, not as distance zero. These selection numbers are
         initial settings for the first pool, not calibrated thresholds.
"""

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "callout_scale_v1" / "batch50"))
from build_panel import SKANI, attrs, components, source_class  # noqa: E402  (same source rules as batch50)

CLONE_ANI, CLONE_AF = 99.99, 95.0
SPECIES_ANI, SPECIES_AF = 95.0, 50.0
N_SPECIES_NEIGHBOURS = 3
AUDIT_MIN, AUDIT_MAX = 10, 20
NO_RUN = dict(host="", host_scientific_name="", isolation_source="", environmental_medium="")


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


def catalog(a):
    data, s2 = a.work / "data", a.work / "s2"
    s2.mkdir(exist_ok=True)
    pool = read_tsv(a.pool_dir / "pool.tsv")
    reports = {}
    for src in ("refseq", "genbank"):
        for line in open(str(a.work / "raw" / ("ncbi_%s.jsonl" % src))):
            d = json.loads(line)
            reports[d["accession"]] = d
    rows, contigs = [], []
    joint = s2 / "joint.fna.gz"
    with gzip.open(str(joint), "wt", compresslevel=6) as out:
        for p in pool:
            acc = p["download_accession"]
            v = json.loads((data / "genomes" / acc / "verify.json").read_text())
            fna = [f["file"] for f in v["files"] if f["file"].endswith("_genomic.fna")][0]
            rep = reports.get(acc) or reports.get(p["genbank_accession"])
            at = attrs(rep)
            cls, note = source_class(at, NO_RUN, rep["assembly_info"]["biosample"].get("package", ""))
            sample = p["biosample"]
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
            rows.append(dict(
                sample_id=sample, assembly=acc, genbank_accession=p["genbank_accession"], role=p["role"],
                organism=p["organism"], genus=p["organism"].split()[0].strip("[]"), tax_id=p["tax_id"], strain=p["strain"],
                bioproject=p["bioproject"], assembly_level=p["assembly_level"], records=v["records"], total_bp=v["fna_bp"],
                plasmid_records=mol.get("Plasmid", 0), annotation=v.get("annotation", ""), linked_read_types=p["linked_read_types"],
                n_usable_runs=p["n_usable_runs"], atypical=p["atypical"], source_class=cls, source_note=note,
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
    write_tsv(a.out / "contig_index.tsv.gz", contigs)
    (s2 / "genomes.list").write_text("".join(str(data / r["fna_path"]) + "\n" for r in rows))
    meta = dict(joint_fna=str(joint), joint_fna_sha256=sha256(joint), samples=len(rows), contigs=len(contigs),
                contig_bp=sum(c["length"] for c in contigs), header="sample_id|assembly|contig",
                check="every joint record SHA256 equals the S1 record SHA256 of its source FNA")
    (a.out / "joint_index.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2), dict(Counter(r["source_class"] for r in rows)))


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


def pairs(a):
    cat = read_tsv(a.out / "catalog.tsv")
    path_to = {str(a.work / "data" / r["fna_path"]): r["sample_id"] for r in cat}
    ids = [r["sample_id"] for r in cat]
    E = {}
    for e in csv.DictReader(open(str(a.work / "s2" / "skani_edges.tsv")), delimiter="\t"):
        x, y = path_to[e["Ref_file"]], path_to[e["Query_file"]]
        E[frozenset((x, y))] = dict(ani=float(e["ANI"]), af={x: float(e["Align_fraction_ref"]), y: float(e["Align_fraction_query"])})
    clone_edges = [tuple(k) for k, v in E.items() if v["ani"] >= CLONE_ANI and min(v["af"].values()) >= CLONE_AF]
    sp_edges = [tuple(k) for k, v in E.items() if v["ani"] >= SPECIES_ANI and min(v["af"].values()) >= SPECIES_AF]
    clone, species = components(ids, clone_edges), components(ids, sp_edges)
    info = {r["sample_id"]: r for r in cat}
    for r in cat:
        r.update(clone_group=clone[r["sample_id"]], species_cluster=species[r["sample_id"]])

    chosen = defaultdict(set)

    def add(x, y, why):
        chosen[frozenset((x, y))].add(why)

    def nbrs(x):
        return sorted(((v["ani"], min(v["af"].values()), next(iter(k - {x}))) for k, v in E.items() if x in k), reverse=True)

    for x in ids:
        nb = nbrs(x)
        cl = [y for _, _, y in nb if clone[y] == clone[x]]
        if cl:
            add(x, cl[0], "clone_partner")
        same = [y for _, _, y in nb if species[y] == species[x] and clone[y] != clone[x]]
        for y in same[:N_SPECIES_NEIGHBOURS]:
            add(x, y, "species_neighbour")
        other = [y for _, _, y in nb if species[y] != species[x]]
        if other:
            add(x, other[0], "nearest_other_species")
    sizes = Counter(species.values())
    linked = Counter(species[r["sample_id"]] for r in cat if r["role"] == "linked")
    audit = sorted((c for c, n in sizes.items() if AUDIT_MIN <= n <= AUDIT_MAX), key=lambda c: (-linked[c], -sizes[c], c))
    audit_cluster = audit[0] if audit else ""
    if audit_cluster:
        for x, y in itertools.combinations(sorted(i for i in ids if species[i] == audit_cluster), 2):
            add(x, y, "audit_all_vs_all")
    rows = []
    for k, why in sorted(chosen.items(), key=lambda kv: sorted(kv[0])):
        x, y = sorted(k)
        e = E.get(k)
        rows.append(dict(
            pair_id="P%05d" % (len(rows) + 1), sample_a=x, sample_b=y, assembly_a=info[x]["assembly"], assembly_b=info[y]["assembly"],
            reasons=";".join(sorted(why)), ani=e["ani"] if e else "", af_a=e["af"][x] if e else "", af_b=e["af"][y] if e else "",
            skani="hit" if e else "no_hit", same_clone=clone[x] == clone[y], same_species_cluster=species[x] == species[y],
            same_genus_name=info[x]["genus"] == info[y]["genus"],
            both_linked=info[x]["role"] == info[y]["role"] == "linked",
            either_linked="linked" in (info[x]["role"], info[y]["role"])))
    write_tsv(a.out / "catalog.tsv", cat)
    write_tsv(a.out / "pairs.tsv", rows)
    groups = [dict(sample_id=r["sample_id"], clone_group=r["clone_group"], species_cluster=r["species_cluster"], genus=r["genus"],
                   role=r["role"]) for r in cat]
    write_tsv(a.out / "groups.tsv", groups)
    no_edge = [i for i in ids if not any(i in k for k in E)]
    summ = dict(
        schema="fna_pool_s2_v1", samples=len(ids), skani_pairs_reported=len(E), possible_pairs=len(ids) * (len(ids) - 1) // 2,
        samples_without_any_skani_hit=no_edge,
        clone_groups=len(set(clone.values())), clone_groups_size_gt1=sum(1 for n in Counter(clone.values()).values() if n > 1),
        species_clusters=len(sizes), species_cluster_sizes=dict(Counter(sizes.values())),
        species_clusters_ge_10=sum(1 for n in sizes.values() if n >= 10), genera=len({r["genus"] for r in cat}),
        source_class=dict(Counter(r["source_class"] for r in cat)),
        source_class_linked=dict(Counter(r["source_class"] for r in cat if r["role"] == "linked")),
        pairs=len(rows), pairs_by_reason=dict(Counter(x for r in rows for x in r["reasons"].split(";"))),
        pairs_both_linked=sum(r["both_linked"] for r in rows), pairs_either_linked=sum(r["either_linked"] for r in rows),
        audit_cluster=dict(id=audit_cluster, size=sizes.get(audit_cluster, 0), linked=linked.get(audit_cluster, 0),
                           organisms=dict(Counter(info[i]["organism"] for i in ids if species[i] == audit_cluster))),
        settings=dict(clone="ANI >= %.2f and AF >= %.0f both" % (CLONE_ANI, CLONE_AF),
                      species="single linkage ANI >= %.0f and AF >= %.0f both" % (SPECIES_ANI, SPECIES_AF),
                      neighbours="1 clone partner, %d nearest same-species non-clone, 1 nearest other-species (skani hit)"
                                 % N_SPECIES_NEIGHBOURS,
                      audit="all pairs of the species cluster with %d-%d members and the most linked samples" % (AUDIT_MIN, AUDIT_MAX),
                      note="initial settings for this pool, not calibrated thresholds; skani reports only AF >= 15"))
    (a.out / "s2_summary.json").write_text(json.dumps(summ, indent=2) + "\n")
    print(json.dumps(summ, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["catalog", "ani", "pairs"])
    p.add_argument("--work", type=Path, required=True, help="pool work dir on scratch (raw/, data/, s2/)")
    p.add_argument("--pool-dir", type=Path, required=True, help="repo dir with pool.tsv")
    p.add_argument("--out", type=Path, required=True, help="repo dir for S2 tables")
    p.add_argument("--threads", type=int, default=16)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    dict(catalog=catalog, pairs=pairs)[a.mode](a) if a.mode != "ani" else ani(a, a.threads)


if __name__ == "__main__":
    main()
