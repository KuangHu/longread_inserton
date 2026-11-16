"""Freeze the metadata-only 50-sample ONT panel (protocol T2 next step) from source/data eligibility.

Inputs come from fetch_metadata.py (ENA runs, exact-BioSample NCBI GenBank reports, downloaded
FASTA stats) and skani. Nothing here reads PILER/RepeatScout output or any candidate/tnp/read count.

Run order (each step caches into --work and says what is missing):
  build_panel.py      -> writes need_gca.txt when FASTA stats are missing
  fetch_metadata.py download --accessions need_gca.txt
  build_panel.py      -> runs skani on the eligible pool, selects, writes need_gcf.txt
  fetch_metadata.py download --accessions need_gcf.txt
  build_panel.py      -> final tables with the GCA/GCF sequence-set check filled in
"""

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SKANI = "/global/scratch/users/kh36969/longread_callout/paraburkholderia_ani/skani"
PANEL_SIZE = 50
CLONE_CAP = 2
CLONE_ANI, CLONE_AF, SPECIES_ANI = 99.99, 95.0, 98.0
MIN_LARGEST, MAX_RECORDS = 3_000_000, 10

ANIMAL_HOST = re.compile(r"homo sapiens|human|patient|sheep|bovine|cattle|mouse|pig\b|horse|chicken|tenualosa|alydus|anasa|"
                         r"jalysus|leptoglossus|riptortus|insect|\bbug\b|mosquito|\btick\b", re.I)
CLINICAL = re.compile(r"sputum|blood|\bpus\b|wound|cystic fibrosis|patient|hospital|clinical|lavage|tracheo|throat|^eta$|"
                      r"urine|abscess|lung|milk|\bgut\b", re.I)
FUNGAL_HOST = re.compile(r"rhizopus", re.I)
PLANT_HOST = re.compile(r"solanum|musa|oryza|arachis|lotus|triticum|zingiber|syzygium|pueraria|rosa\b|boxwood|tomato|rice|"
                        r"turmeric|ginger|mimosa|zea\b|maize|glycine|phaseolus|medicago|arabidopsis|citrus|vitis|nicotiana", re.I)
ENRICH = re.compile(r"enrich|consortium", re.I)
BUILT = re.compile(r"industrial|pharmaceutical|\biss\b|spacecraft|clean ?room", re.I)
PLANT_TEXT = re.compile(r"root|rhizosph|nodule|\bleaf\b|leaves|xylem|\bplant\b|seed|stem|potato|tomato|banana|rice|wheat|maize|"
                        r"citri|vetiver|juniperus|\barum\b|fruit|ginger|turmeric", re.I)
ENV_TEXT = re.compile(r"soil|sediment|water|sludge|\bmine\b|tar pit|yard|forest|lake|river|pond|groundwater|environmental|"
                      r"field|hyphosphere|ENVO", re.I)
MISSING = re.compile(r"^(|missing|not collected|not provided|not applicable|na|n/a|unknown)$", re.I)
PCR = re.compile(r"PCR", re.I)
ALREADY_RUN = {"acceptance": HERE.parent / "test_panel.tsv", "batch10": HERE.parent / "batch10" / "samples.tsv"}


def attrs(report):
    return {a["name"]: a.get("value", "") for a in report["assembly_info"]["biosample"].get("attributes", [])}


def clean(v):
    return "" if MISSING.match((v or "").strip()) else v.strip()


def source_class(a, run, package):
    """Fixed keyword rules over BioSample host / isolation_source / env_medium (ENA as fallback)."""
    host = clean(a.get("host")) or clean(run["host"]) or clean(run["host_scientific_name"])
    iso = clean(a.get("isolation_source")) or clean(run["isolation_source"])
    med = clean(a.get("env_medium")) or clean(run["environmental_medium"])
    stype = clean(a.get("sample_type"))
    text = " ".join([iso, med])
    if host and ANIMAL_HOST.search(host):
        return "clinical/animal", "animal or human host: " + host
    if CLINICAL.search(text) or CLINICAL.search(host):
        return "clinical/animal", "clinical text: " + (text.strip() or host)
    if re.search(r"mixed culture", stype, re.I):
        return "unclear", "sample_type mixed culture (not a pure isolate)"
    if host and FUNGAL_HOST.search(host):
        return "unclear", "fungal endosymbiont; host " + host + (" / " + iso if iso else "")
    if host and PLANT_HOST.search(host):
        return "plant", "plant host: " + host + (" / " + text.strip() if text.strip() else "")
    if ENRICH.search(text):
        return "enrichment", iso or med
    if BUILT.search(text):
        return "unclear", "built/industrial: " + text.strip()
    if PLANT_TEXT.search(text):
        return "plant", text.strip()
    if ENV_TEXT.search(text):
        return "environmental", text.strip()
    if host:
        return "unclear", "unrecognised host: " + host
    if package.startswith("Pathogen.cl"):
        return "clinical/animal", "Pathogen.cl package with no host/source"
    return "unclear", "no host/source metadata"


def library_class(run):
    if run["library_strategy"] != "WGS":
        return "excluded_library_strategy"
    sel = run["library_selection"]
    if PCR.search(sel):
        return "pcr"
    if sel in ("RANDOM", "other", "unspecified", "size fractionation"):
        return "non_pcr"
    return "excluded_library_selection"


def load_tsv(path):
    with open(path) as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as h:
        w = csv.DictWriter(h, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def components(nodes, edges):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        parent[find(a)] = find(b)
    groups = defaultdict(list)
    for n in nodes:
        groups[find(n)].append(n)
    return {n: sorted(g)[0] for g in groups.values() for n in g}


def sample_id(report, run):
    names = report["organism"].get("infraspecific_names", {})
    name = names.get("strain") or names.get("isolate") or run["strain"] or run["isolate"] \
        or report["assembly_info"]["biosample"]["accession"]
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout/batch50_meta"))
    p.add_argument("--out", type=Path, default=HERE)
    a = p.parse_args()
    raw, gen = a.work / "raw", a.work / "genomes"

    runs = load_tsv(raw / "ena_ont_runs.tsv")
    reports = [json.loads(l) for l in open(raw / "ncbi_biosample_reports.jsonl")]
    gca = defaultdict(list)
    for r in reports:
        if r["source_database"] == "SOURCE_DATABASE_GENBANK" and r["assembly_info"]["assembly_status"] == "current":
            gca[r["assembly_info"]["biosample"]["accession"]].append(r)
    done = {}
    for ref in ALREADY_RUN.values():
        for r in load_tsv(ref):
            done[r["run"]] = done[r["biosample"]] = done[r.get("assembly") or r.get("gca")] = r["sample_id"]

    genomic = [r for r in runs if r["library_source"] == "GENOMIC"]
    rows = []
    for run in genomic:
        matches = gca.get(run["sample_accession"], [])
        if not matches:
            rows.append(dict(run=run["run_accession"], biosample=run["sample_accession"], organism=run["scientific_name"],
                             library_strategy=run["library_strategy"], library_selection=run["library_selection"],
                             status="excluded", reason="no_exact_biosample_genbank_assembly"))
            continue
        for rep in matches:
            ai, st = rep["assembly_info"], rep["assembly_stats"]
            at = attrs(rep)
            package = ai["biosample"].get("package", "")
            cls, note = source_class(at, run, package)
            org = rep["organism"]["organism_name"]
            rows.append(dict(
                run=run["run_accession"], biosample=run["sample_accession"], study=run["study_accession"],
                bioproject=ai.get("bioproject_accession", ""), gca=rep["accession"],
                gcf=rep.get("paired_accession", ""), organism=org, genus=org.split()[0],
                strain=rep["organism"].get("infraspecific_names", {}).get("strain") or
                rep["organism"].get("infraspecific_names", {}).get("isolate") or run["strain"] or run["isolate"],
                sample_id=sample_id(rep, run), assembly_level=ai["assembly_level"],
                report_records=st.get("number_of_scaffolds", ""), report_total_bp=st.get("total_sequence_length", ""),
                instrument_model=run["instrument_model"], library_strategy=run["library_strategy"],
                library_source=run["library_source"], library_selection=run["library_selection"],
                read_count=run["read_count"], base_count=run["base_count"], fastq_bytes=run["fastq_bytes"],
                fastq_md5=run["fastq_md5"], fastq_ftp=run["fastq_ftp"], biosample_package=package,
                host=clean(at.get("host")) or clean(run["host"]),
                isolation_source=clean(at.get("isolation_source")) or clean(run["isolation_source"]),
                env_medium=clean(at.get("env_medium")) or clean(run["environmental_medium"]),
                geo_loc_name=clean(at.get("geo_loc_name")) or clean(run["country"]),
                collection_date=clean(at.get("collection_date")) or clean(run["collection_date"]),
                source_class=cls, source_note=note, library_class=library_class(run),
                genus_tier="1_paraburkholderia" if org.startswith("Paraburkholderia") else "2_burkholderiaceae",
            ))

    # Assembly continuity is judged from the downloaded GenBank FASTA, with the report only as a prefilter.
    def prefilter(r):
        return r["assembly_level"] == "Complete Genome" or int(r["report_records"] or 99) <= MAX_RECORDS
    excluded_species = lambda r: re.match(r"Burkholderia (pseudo)?mallei\b", r["organism"])
    need = sorted({r["gca"] for r in rows if "gca" in r and prefilter(r) and not excluded_species(r)}
                  | {k for k in done if k.startswith("GCA_")})
    missing = [g for g in need if not (gen / (g + ".json")).exists()]
    if missing:
        (a.work / "need_gca.txt").write_text("\n".join(missing) + "\n")
        raise SystemExit("%d GCA FASTA missing; run fetch_metadata.py download --accessions %s"
                         % (len(missing), a.work / "need_gca.txt"))
    fa = {g: json.loads((gen / (g + ".json")).read_text()) for g in need}

    for r in rows:
        if r.get("status"):
            continue
        f = fa.get(r["gca"])
        if f:
            r.update(records=f["records"], total_bp=f["total_bp"], largest_bp=f["largest_bp"], n50_bp=f["n50_bp"],
                     n_bases=f["n_bases"], gca_fasta_sha256=f["fasta_sha256"],
                     bases_per_assembly_bp="%.1f" % (int(r["base_count"] or 0) / f["total_bp"]))
        complete = r["assembly_level"] == "Complete Genome"
        contiguous = f and (complete or (f["records"] <= MAX_RECORDS and f["largest_bp"] >= MIN_LARGEST))
        r["continuity"] = "complete" if complete and f else "high_contiguity_draft" if contiguous else "fail"
        if r["library_strategy"] != "WGS":
            r.update(status="excluded", reason="library_strategy_not_WGS")
        elif excluded_species(r):
            r.update(status="excluded", reason="excluded_species_B_pseudomallei_or_mallei")
        elif not (r["fastq_md5"] and r["fastq_bytes"]) or ";" in r["fastq_ftp"]:
            r.update(status="excluded", reason="ena_fastq_bytes_md5_missing_or_not_single_file")
        elif not contiguous:
            r.update(status="excluded", reason="continuity_fail" + (
                ": %s records, largest %s bp" % (f["records"], f["largest_bp"]) if f else
                ": %s, %s records in report" % (r["assembly_level"], r["report_records"])))
        elif r["run"] in done or r["biosample"] in done or r["gca"] in done:
            r.update(status="excluded", reason="already_run: " + done.get(r["run"]) if r["run"] in done
                     else "already_run: " + done.get(r["biosample"], done.get(r["gca"], "")))
        elif r["source_class"] == "clinical/animal":
            r.update(status="excluded", reason="source_clinical_or_animal: " + r["source_note"])
        elif r["library_class"] == "excluded_library_selection":
            r.update(status="excluded", reason="library_selection_not_genomic_random_like: " + r["library_selection"])
        elif r["source_note"].startswith("sample_type mixed culture"):
            r.update(status="reserve", reason="not_pure_isolate: " + r["source_note"])
        elif r["library_class"] == "pcr":
            r.update(status="pcr_reserve", reason="PCR library (" + r["library_selection"] + "): stress/reserve only")
        else:
            r.update(status="eligible", reason="")

    # One run per BioSample: non-PCR, RANDOM, then more bases (data volume only).
    pool = [r for r in rows if r.get("status") in ("eligible", "pcr_reserve")]
    by_bs = defaultdict(list)
    for r in pool:
        by_bs[r["biosample"]].append(r)
    for bs, rs in by_bs.items():
        rs.sort(key=lambda r: (r["status"] != "eligible", r["library_selection"] != "RANDOM", -int(r["base_count"] or 0), r["run"]))
        for r in rs[1:]:
            r.update(status="reserve", reason="alternate_run_same_biosample: chosen " + rs[0]["run"])
    # A run whose BioSample has two current GCAs keeps the more contiguous one.
    by_gca_bs = defaultdict(list)
    for r in rows:
        if r.get("status") in ("eligible", "pcr_reserve"):
            by_gca_bs[(r["run"], r["biosample"])].append(r)
    for rs in by_gca_bs.values():
        rs.sort(key=lambda r: (r["continuity"] != "complete", r["records"], r["gca"]))
        for r in rs[1:]:
            r.update(status="reserve", reason="second_current_gca_same_biosample: chosen " + rs[0]["gca"])

    # skani over the remaining pool plus the 15 already-run assemblies.
    pool = [r for r in rows if r.get("status") in ("eligible", "pcr_reserve")]
    ani_set = sorted({r["gca"] for r in pool} | {k for k in done if k.startswith("GCA_")})
    ani_dir = a.work / "ani"
    ani_dir.mkdir(exist_ok=True)
    lst = ani_dir / "genomes.txt"
    lst_text = "\n".join(str(gen / (g + ".fna")) for g in ani_set) + "\n"
    key = hashlib.sha256(lst_text.encode()).hexdigest()[:16]
    edges_path = ani_dir / ("triangle_%s.tsv" % key)
    cmd = [SKANI, "triangle", "-l", str(lst), "-E", "--slow", "-t", "8", "-o", str(edges_path)]
    if not edges_path.exists():
        lst.write_text(lst_text)
        subprocess.run(cmd, check=True)
    acc = lambda path: Path(path).name.replace(".fna", "")
    edges = load_tsv(edges_path)
    clone_e, species_e, ani95_e = [], [], []
    for e in edges:
        q, rf, ani = acc(e["Query_file"]), acc(e["Ref_file"]), float(e["ANI"])
        if ani >= CLONE_ANI and min(float(e["Align_fraction_ref"]), float(e["Align_fraction_query"])) >= CLONE_AF:
            clone_e.append((q, rf))
        if ani >= SPECIES_ANI:
            species_e.append((q, rf))
        if ani >= 95:
            ani95_e.append((q, rf))
    clone = components(ani_set, clone_e)
    species = components(ani_set, species_e)
    ani95 = components(ani_set, ani95_e)
    done_gca = {k: v for k, v in done.items() if k.startswith("GCA_")}
    for r in rows:
        if r.get("gca") in clone:
            r["clone_group"] = "C_" + clone[r["gca"]]
            r["species_cluster"] = "S_" + species[r["gca"]]
            r["ani95_group"] = "G_" + ani95[r["gca"]]
            same = sorted(v for g, v in done_gca.items() if clone[g] == clone[r["gca"]])
            r["clone_already_run"] = ",".join(same)

    # Rank before any discovery: genus tier, RANDOM, source class, Complete, read bases per assembly bp.
    src_rank = {"environmental": 0, "plant": 0, "enrichment": 0, "unclear": 1}
    keyf = lambda r: (r["genus_tier"], r["library_selection"] != "RANDOM", src_rank[r["source_class"]],
                      r["continuity"] != "complete", -float(r["bases_per_assembly_bp"]), r["run"])
    eligible = sorted([r for r in rows if r.get("status") == "eligible"], key=keyf)
    for i, r in enumerate(eligible, 1):
        r["rank"] = i
    # Diversity after clone collapse: pass 1 one per species cluster, pass 2 one per clone group,
    # pass 3 second clone member; Paraburkholderia rows go before family fill in every pass.
    # Already-run members count against the clone cap and as already-represented species clusters.
    used = Counter(clone[g] for g in done_gca)
    seen_species = {species[g] for g in done_gca}
    selected = []
    for pass_no in (1, 2, 3):
        for tier in ("1_paraburkholderia", "2_burkholderiaceae"):
            for r in eligible:
                if len(selected) >= PANEL_SIZE or r.get("selected") or r["genus_tier"] != tier:
                    continue
                c, s = clone[r["gca"]], species[r["gca"]]
                if pass_no == 1 and (s in seen_species or used[c]):
                    continue
                if pass_no == 2 and used[c]:
                    continue
                if pass_no == 3 and used[c] >= CLONE_CAP:
                    continue
                used[c] += 1
                seen_species.add(s)
                r.update(selected=True, selection_pass=pass_no)
                selected.append(r)
    for r in eligible:
        if r.get("selected"):
            r.update(status="selected", reason="pass %d (%s)" % (r["selection_pass"], {
                1: "first member of a species cluster", 2: "first member of a clone group",
                3: "second member of a clone group"}[r["selection_pass"]]))
        elif used[clone[r["gca"]]]:
            # Near-clone of a selected or already-run sample: kept as clone reserve, not deleted.
            r.update(status="clone_reserve", reason="clone group %s has %d selected/already-run (cap %d; already run: %s)" % (
                r["clone_group"], used[clone[r["gca"]]], CLONE_CAP, r["clone_already_run"] or "none"))
        else:
            r.update(status="ranked_reserve", reason="independent clone group beyond panel size")
    # Order inside a clone group, for the uniqueness check: already-run members count first.
    order = Counter(clone[g] for g in done_gca)
    for r in sorted(selected, key=lambda r: r["rank"]):
        order[clone[r["gca"]]] += 1
        r["clone_member"] = "first" if order[clone[r["gca"]]] == 1 else "second_allowed"

    # GCA/GCF sequence-set check for the selected rows and the first 20 independent reserves.
    check_rows = selected + [r for r in eligible if r["status"] == "ranked_reserve"][:20]
    need_gcf = sorted({r["gcf"] for r in check_rows if r["gcf"] and not (gen / (r["gcf"] + ".json")).exists()})
    (a.work / "need_gcf.txt").write_text("\n".join(need_gcf) + ("\n" if need_gcf else ""))
    for r in rows:
        if not r.get("gca") or r.get("status") not in ("selected", "clone_reserve", "ranked_reserve", "pcr_reserve"):
            continue
        if not r["gcf"]:
            r["gca_gcf_check"] = "no_gcf_use_gca"
        elif (gen / (r["gcf"] + ".json")).exists():
            g = json.loads((gen / (r["gcf"] + ".json")).read_text())
            r["gca_gcf_check"] = "done_identical" if g["record_sha256"] == fa[r["gca"]]["record_sha256"] else "done_different_use_gca"
        else:
            r["gca_gcf_check"] = "pending"

    write_outputs(a, rows, genomic, runs, reports, selected, eligible, edges_path, cmd, ani_set, clone, species, done_gca, ani95)


FUNNEL = [("library_strategy_not_WGS", "WGS"), ("excluded_species", "not B. pseudomallei/mallei"),
          ("ena_fastq", "ENA single FASTQ with bytes and MD5"), ("continuity_fail", "continuity"),
          ("already_run", "not already run"), ("source_clinical_or_animal", "source not clinical/animal"),
          ("library_selection_not", "library selection genomic-random-like"), ("not_pure_isolate", "pure isolate"),
          ("PCR library", "non-PCR (PCR -> pcr_reserve)"), ("alternate_run_same_biosample", "one run per BioSample"),
          ("second_current_gca", "one assembly per BioSample")]


def funnel(matched):
    """Sequential counts of exact-BioSample run/assembly rows, in the order the filters are applied."""
    left, out = len(matched), [dict(step="exact BioSample run x current GCA", rows=len(matched))]
    for prefix, label in FUNNEL:
        n = sum(1 for r in matched if r["reason"].startswith(prefix))
        left -= n
        out.append(dict(step=label, removed=n, rows=left))
    return out


COLS = ["rank", "sample_id", "organism", "strain", "genus_tier", "run", "biosample", "bioproject", "study", "gca", "gcf",
        "gca_gcf_check", "assembly_level", "continuity", "records", "total_bp", "largest_bp", "n50_bp", "n_bases",
        "instrument_model", "library_strategy", "library_source", "library_selection", "read_count", "base_count",
        "bases_per_assembly_bp", "fastq_bytes", "fastq_md5", "fastq_ftp", "host", "isolation_source", "env_medium",
        "geo_loc_name", "collection_date", "biosample_package", "source_class", "source_note", "clone_group",
        "clone_member", "clone_already_run", "species_cluster", "ani95_group", "gca_fasta_sha256", "status", "reason"]


def write_outputs(a, rows, genomic, runs, reports, selected, eligible, edges_path, cmd, ani_set, clone, species, done_gca, ani95):
    out = a.out
    sel = sorted(selected, key=lambda r: r["rank"])
    for i, r in enumerate(sel, 1):
        r["order"] = i
    write_tsv(out / "selection_manifest.tsv", sel, ["order"] + COLS)
    write_tsv(out / "ranked_candidates.tsv", eligible, COLS + ["selection_pass"])
    write_tsv(out / "clone_reserve.tsv", [r for r in eligible if r["status"] == "clone_reserve"], COLS)
    excl_cols = ["run", "biosample", "gca", "organism", "assembly_level", "records", "largest_bp", "library_strategy",
                 "library_selection", "source_class", "source_note", "status", "reason"]
    write_tsv(out / "exclusion_reasons.tsv", sorted([r for r in rows if r["status"] not in ("selected",) and r not in eligible],
                                                    key=lambda r: (r["run"], r.get("gca", ""))), excl_cols)
    groups = []
    members = defaultdict(list)
    for g in ani_set:
        members[clone[g]].append(g)
    info = {r["gca"]: r for r in rows if r.get("gca")}
    for c, gs in sorted(members.items()):
        sel_n = sum(1 for g in gs if info.get(g, {}).get("status") == "selected")
        for g in sorted(gs):
            r = info.get(g, {})
            groups.append(dict(clone_group="C_" + c, species_cluster="S_" + species[g], ani95_group="G_" + ani95[g], gca=g, organism=r.get("organism", ""),
                               sample_id=r.get("sample_id", done_gca.get(g, "")), run=r.get("run", ""),
                               status="already_run" if g in done_gca else r.get("status", ""),
                               group_size=len(gs), group_selected=sel_n,
                               group_already_run=sum(1 for x in gs if x in done_gca)))
    write_tsv(out / "clone_groups.tsv", groups, list(groups[0]))

    stage = Counter()
    gen_runs = len(genomic)
    matched = [r for r in rows if r.get("gca")]
    summary = dict(
        schema="batch50_selection_v1", date=json.loads((a.work / "raw/query_provenance.json").read_text())["date"],
        ena_ont_runs_all_sources=len(runs), ena_ont_genomic_runs_scanned=gen_runs,
        ena_ont_genomic_biosamples=len({r["sample_accession"] for r in genomic}),
        exact_biosample_run_assembly_rows=len(matched),
        exact_biosample_runs=len({r["run"] for r in matched}),
        exact_biosample_assemblies=len({r["gca"] for r in matched}),
        funnel=funnel(matched),
        status_counts=dict(Counter(r["status"] for r in rows)),
        exclusion_reason_counts=dict(Counter(r["reason"].split(":")[0] for r in rows if r["status"] in ("excluded", "reserve"))),
        eligible_new_samples=len(eligible),
        eligible_by_genus_tier=dict(Counter(r["genus_tier"] for r in eligible)),
        eligible_clone_groups=len({clone[r["gca"]] for r in eligible}),
        eligible_species_clusters=len({species[r["gca"]] for r in eligible}),
        eligible_clone_groups_not_already_run=len({clone[r["gca"]] for r in eligible if not r["clone_already_run"]}),
        max_new_samples_under_clone_cap=None,
        selected=len(selected), selected_by_pass=dict(Counter(r["selection_pass"] for r in selected)),
        selected_clone_groups=len({clone[r["gca"]] for r in selected}),
        selected_species_clusters=len({species[r["gca"]] for r in selected}),
        selected_ani95_groups=len({r["ani95_group"] for r in selected}),
        selected_ani95_groups_shared_with_already_run=len({r["ani95_group"] for r in selected} & {"G_" + ani95[g] for g in done_gca}),
        selected_organism_names=len({" ".join(r["organism"].split()[:2]) for r in selected}),
        eligible_ani95_groups=len({r["ani95_group"] for r in eligible}),
        selected_by_genus=dict(Counter(r["genus"] for r in selected)),
        selected_by_source=dict(Counter(r["source_class"] for r in selected)),
        selected_by_continuity=dict(Counter(r["continuity"] for r in selected)),
        selected_by_library_selection=dict(Counter(r["library_selection"] for r in selected)),
        selected_gca_gcf_check=dict(Counter(r.get("gca_gcf_check", "") for r in selected)),
        selected_fastq_bytes_total=sum(int(r["fastq_bytes"]) for r in selected),
        selected_base_count_total=sum(int(r["base_count"]) for r in selected),
        panel_size_target=PANEL_SIZE, fifty_available=len(selected) >= PANEL_SIZE,
    )
    used = Counter(clone[g] for g in done_gca)
    cap_room = Counter()
    for r in eligible:
        c = clone[r["gca"]]
        if used[c] + cap_room[c] < CLONE_CAP:
            cap_room[c] += 1
    summary["max_new_samples_under_clone_cap"] = sum(cap_room.values())
    (out / "selection_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    prov = json.loads((a.work / "raw/query_provenance.json").read_text())
    prov.update(
        skani=dict(command=cmd, version=subprocess.run([SKANI, "--version"], capture_output=True, text=True).stdout.strip(),
                   binary_sha256=hashlib.sha256(Path(SKANI).read_bytes()).hexdigest(), genomes=len(ani_set),
                   edges_file=str(edges_path), edges_sha256=hashlib.sha256(edges_path.read_bytes()).hexdigest()),
        rules=dict(clone="single linkage ANI >= %.2f and min(AF_ref, AF_query) >= %.0f" % (CLONE_ANI, CLONE_AF),
                   species_cluster="single linkage ANI >= %.0f (skani reports pairs with AF >= 15)" % SPECIES_ANI,
                   ani95_group="single linkage ANI >= 95 (species-level reference only; not used for selection)",
                   clone_cap="%d per clone group in the primary panel; already-run members count against the cap" % CLONE_CAP,
                   continuity="Complete Genome, or GenBank FASTA records <= %d and largest record >= %d bp" % (MAX_RECORDS, MIN_LARGEST),
                   library="WGS/GENOMIC; RANDOM, other, unspecified, size fractionation = non-PCR; any *PCR* = pcr_reserve; others excluded",
                   species_excluded="Burkholderia pseudomallei, Burkholderia mallei",
                   source="fixed keyword rules in build_panel.source_class over BioSample host/isolation_source/env_medium; ENA fallback",
                   ranking="genus tier (Paraburkholderia first), RANDOM first, environmental/plant/enrichment before unclear, "
                           "Complete before draft, read bases per assembly bp descending, run accession",
                   diversity="pass 1 first member per ANI>=98 cluster not already run, pass 2 first member per clone group, pass 3 second clone member",
                   not_used="no candidate, tnp, read-spanned or other PILER/RepeatScout output"),
        raw_inputs={f: hashlib.sha256((a.work / "raw" / f).read_bytes()).hexdigest()
                    for f in ("ena_ont_runs.tsv", "ncbi_biosample_reports.jsonl")},
        raw_copies_in_repo="raw/*.gz (gzip of the scratch files above and of the skani edge list)",
        downloads="GenBank (and paired RefSeq for checked rows) genomic FASTA via NCBI Datasets v2 download API; per-record SHA256",
        script_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE / "build_panel.py", HERE / "fetch_metadata.py")},
    )
    (out / "query_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "exclusion_reason_counts"}, indent=2))
    print("reasons", json.dumps(summary["exclusion_reason_counts"], indent=1))


if __name__ == "__main__":
    main()
