"""S0: assembly catalog, ONT/PacBio run list, assembly-run links and the S1 download manifest for one taxon.

Metadata only (docs/PLAN_CROSS_ASSEMBLY_GFF3_SHORT_INSERTIONS.md, S0). No reads and no FASTA are downloaded.

fetch  NCBI Datasets genome summaries (RefSeq and GenBank, current assemblies) for --taxon; ENA read_run
       rows for the taxon on OXFORD_NANOPORE and PACBIO_SMRT; SRA XML of the PacBio runs (original submitted
       file names, used as read-type evidence). Raw responses go to <work>/raw with SHA256 in provenance.json.
build  Joins run sample_accession to assembly BioSample exactly (link_basis biosample_matched; source-run
       evidence is not checked here). One assembly per BioSample: the GenBank assembly, its RefSeq pair as the
       download target when NCBI does not flag them as different. Pool = every BioSample with a usable
       genomic run (role linked) + every other RefSeq Complete Genome (role context_refseq_complete).
       Writes the tables below into --out. GCA/GCF identity is only checked after download (S1).

Python 3.6+, standard library only (plus the datasets CLI and curl for fetch).
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import time
import urllib.parse
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

DATASETS = "/global/home/users/kh36969/.conda/envs/claude-env/bin/datasets"
ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RUN_FIELDS = ("run_accession,experiment_accession,sample_accession,secondary_sample_accession,study_accession,"
              "tax_id,scientific_name,instrument_platform,instrument_model,library_strategy,library_source,"
              "library_selection,library_layout,library_name,read_count,base_count,fastq_ftp,fastq_bytes,fastq_md5,"
              "submitted_ftp,submitted_bytes,submitted_md5,submitted_format,first_public")
PLATFORMS = ("OXFORD_NANOPORE", "PACBIO_SMRT")
LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
HIFI = re.compile(r"\bccs\b|hifi|\.q20\b|reads\.(?:fastq|fq|bam)|\.ccs\.", re.I)
MIN_MEAN_LENGTH = 500  # shorter mean = short-read data under a long-read label (seen: 150-302 bp R1/R2 files)
SUBREADS = re.compile(r"subreads?\.bam|\.subreads\.|\.bax\.h5|\.bas\.h5", re.I)


def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def curl(url, out):
    subprocess.run(["curl", "-sS", "--fail", "--retry", "5", "--retry-delay", "10", "-o", str(out), url], check=True)


def fetch(a):
    raw = a.work / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    prov = dict(date=time.strftime("%Y-%m-%dT%H:%M:%S%z"), taxon=a.taxon, files={})
    for src in ("refseq", "genbank"):
        out = raw / ("ncbi_%s.jsonl" % src)
        cmd = [DATASETS, "summary", "genome", "taxon", str(a.taxon), "--assembly-source", src, "--as-json-lines"]
        with open(str(out), "w") as h:
            subprocess.run(cmd, stdout=h, stderr=subprocess.DEVNULL, check=True)
        prov["files"][out.name] = dict(command=" ".join(cmd[1:]), rows=sum(1 for _ in open(str(out))), sha256=sha256(out))
    q = "tax_tree(%d) AND (%s)" % (a.taxon, " OR ".join('instrument_platform="%s"' % p for p in PLATFORMS))
    url = ENA + "?" + urllib.parse.urlencode(dict(result="read_run", query=q, fields=RUN_FIELDS, format="tsv", limit=0))
    out = raw / "ena_runs.tsv"
    curl(url, out)
    runs = list(csv.DictReader(open(str(out)), delimiter="\t"))
    prov["files"][out.name] = dict(url=url, rows=len(runs), sha256=sha256(out))
    pb = sorted(r["run_accession"] for r in runs if r["instrument_platform"] == "PACBIO_SMRT")
    for i in range(0, len(pb), 100):
        out = raw / ("sra_pacbio_%03d.xml" % (i // 100))
        url = EUTILS + "?" + urllib.parse.urlencode(dict(db="sra", id=",".join(pb[i:i + 100]), rettype="xml"))
        curl(url, out)
        prov["files"][out.name] = dict(url=url, runs=len(pb[i:i + 100]), sha256=sha256(out))
        time.sleep(1)
    (raw / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps({k: v.get("rows", v.get("runs")) for k, v in prov["files"].items()}, indent=2))


def sra_files(raw):
    """run -> original submitted file names from SRA XML (SRAFile entries that are not SRA-made)."""
    out = defaultdict(list)
    for p in sorted(raw.glob("sra_pacbio_*.xml")):
        for run in ET.parse(str(p)).getroot().iter("RUN"):
            acc = run.get("accession")
            for f in run.iter("SRAFile"):
                name = f.get("filename") or ""
                if name and name != acc and not name.startswith(acc + "."):
                    out[acc].append(name)
    return out


def read_type(r, submitted):
    """Read type and the evidence used; ENA FASTQ names (SRA-made, often *_subreads) are not evidence."""
    if r["instrument_platform"] == "OXFORD_NANOPORE":
        return "ont", "instrument_platform"
    names = submitted.get(r["run_accession"], []) + [x.rsplit("/", 1)[-1] for x in r["submitted_ftp"].split(";") if x]
    text = " ".join(names)
    if HIFI.search(text) and not SUBREADS.search(text):
        return "pacbio_hifi_ccs", "submitted_file:" + text[:200]
    if SUBREADS.search(text) and not HIFI.search(text):
        return "pacbio_clr_subreads", "submitted_file:" + text[:200]
    return "pacbio_unresolved", ("submitted_file:" + text[:200]) if text else "no_submitted_file_name"


def usable(r):
    """Genomic WGS, non-PCR, long reads (mean >= MIN_MEAN_LENGTH), with an ENA FASTQ that has MD5 and bytes."""
    reasons = []
    if r["library_source"] != "GENOMIC":
        reasons.append("library_source=" + r["library_source"])
    if r["library_strategy"] != "WGS":
        reasons.append("library_strategy=" + r["library_strategy"])
    if "PCR" in r["library_selection"].upper():
        reasons.append("library_selection=" + r["library_selection"])
    if not (r["fastq_md5"] and r["fastq_bytes"]):
        reasons.append("no_ena_fastq_md5")
    if r["instrument_model"] == "Onso":
        reasons.append("short_read_instrument=Onso")
    if r["read_count"] and mean_length(r) < MIN_MEAN_LENGTH:
        reasons.append("mean_read_length<%d" % MIN_MEAN_LENGTH)
    return reasons


def mean_length(r):
    n, b = int(r["read_count"] or 0), int(r["base_count"] or 0)
    return b / n if n else 0.0


def load_assemblies(raw):
    rows = {}
    for src in ("refseq", "genbank"):
        for line in open(str(raw / ("ncbi_%s.jsonl" % src))):
            d = json.loads(line)
            ai, st = d["assembly_info"], d.get("assembly_stats", {})
            pair = ai.get("paired_assembly", {})
            ann = d.get("annotation_info", {})
            rows[d["accession"]] = dict(
                accession=d["accession"], source=src, paired_accession=pair.get("accession", ""),
                pair_status=pair.get("status", ""), ncbi_flags_pair_different=bool(pair.get("refseq_genbank_are_different")),
                biosample=ai.get("biosample", {}).get("accession", ""), bioproject=ai.get("bioproject_accession", ""),
                organism=d.get("organism", {}).get("organism_name", ""), tax_id=d.get("organism", {}).get("tax_id", ""),
                strain=d.get("organism", {}).get("infraspecific_names", {}).get("strain", ""),
                assembly_level=ai.get("assembly_level", ""), assembly_status=ai.get("assembly_status", ""),
                release_date=ai.get("release_date", ""), sequencing_tech=ai.get("sequencing_tech", ""),
                assembly_method=ai.get("assembly_method", ""),
                atypical=";".join(ai.get("atypical", {}).get("warnings", [])),
                total_bp=st.get("total_sequence_length", ""), records=st.get("number_of_component_sequences", ""),
                annotation_provider=ann.get("provider", ""), annotation_name=ann.get("name", ""),
                annotation_date=ann.get("release_date", ""),
                protein_coding=ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", ""),
                checkm_completeness=d.get("checkm_info", {}).get("completeness", ""),
                checkm_contamination=d.get("checkm_info", {}).get("contamination", ""))
    return rows


def choose(gca_rows, asm):
    """One GenBank assembly per BioSample: current, then level, then newest release, then accession."""
    cur = [r for r in gca_rows if r["assembly_status"] == "current"] or gca_rows
    best = sorted(cur, key=lambda r: (LEVEL_RANK.get(r["assembly_level"], 9), "".join(chr(255 - ord(c)) for c in r["release_date"]),
                                      r["accession"]))[0]
    gcf = asm.get(best["paired_accession"])
    use_gcf = bool(gcf) and gcf["source"] == "refseq" and not best["ncbi_flags_pair_different"]
    return best, (gcf if use_gcf else None)


def write_tsv(path, rows, cols=None):
    cols = cols or list(rows[0])
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "wt") as h:
        w = csv.DictWriter(h, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def build(a):
    raw = a.work / "raw"
    a.out.mkdir(parents=True, exist_ok=True)
    asm = load_assemblies(raw)
    runs = list(csv.DictReader(open(str(raw / "ena_runs.tsv")), delimiter="\t"))
    submitted = sra_files(raw)
    by_bs = defaultdict(list)
    for r in asm.values():
        if r["source"] == "genbank" and r["biosample"]:
            by_bs[r["biosample"]].append(r)
    # RefSeq rows whose GenBank pair is missing from the GenBank listing still count for their BioSample.
    for r in asm.values():
        if r["source"] == "refseq" and r["biosample"] and r["paired_accession"] not in asm:
            by_bs[r["biosample"]].append(r)

    run_rows, links = [], []
    for r in runs:
        rt, ev = read_type(r, submitted)
        why = usable(r)
        cands = by_bs.get(r["sample_accession"], [])
        run_rows.append(dict(r, read_type=rt, read_type_evidence=ev, mean_read_length=round(mean_length(r)), usable=not why, unusable_reasons=";".join(why),
                             submitted_names_sra=";".join(submitted.get(r["run_accession"], [])),
                             n_assemblies_same_biosample=len(cands)))
        for c in cands:
            links.append(dict(assembly_accession=c["accession"], paired_accession=c["paired_accession"],
                              biosample=r["sample_accession"], run=r["run_accession"], study=r["study_accession"],
                              read_type=rt, link_basis="biosample_matched", source_run_evidence="not_checked",
                              run_usable=not why, run_tax_id=r["tax_id"], run_organism=r["scientific_name"],
                              assembly_organism=c["organism"],
                              genus_agrees=r["scientific_name"].split()[:1] == c["organism"].split()[:1],
                              assembly_sequencing_tech=c["sequencing_tech"]))
    linked_bs = defaultdict(set)
    for l in links:
        if l["run_usable"]:
            linked_bs[l["biosample"]].add(l["read_type"])

    pool, seen = [], set()
    for bs in sorted(linked_bs):
        best, gcf = choose(by_bs[bs], asm)
        pool.append((bs, best, gcf, "linked"))
        seen.add(bs)
    for r in sorted(asm.values(), key=lambda r: r["accession"]):
        if r["source"] == "refseq" and r["assembly_level"] == "Complete Genome" and r["biosample"] not in seen \
                and r["assembly_status"] == "current":
            best, gcf = choose(by_bs.get(r["biosample"]) or [r], asm)
            pool.append((r["biosample"], best, gcf, "context_refseq_complete"))
            seen.add(r["biosample"])

    pool_rows, manifest = [], []
    for bs, best, gcf, role in pool:
        target = gcf or best
        alts = sorted(x["accession"] for x in by_bs.get(bs, []) if x["accession"] != best["accession"])
        rtypes = sorted(linked_bs.get(bs, []))
        pool_rows.append(dict(
            biosample=bs, role=role, download_accession=target["accession"],
            genbank_accession=best["accession"] if best["source"] == "genbank" else best["paired_accession"],
            refseq_accession=gcf["accession"] if gcf else ("" if best["source"] == "genbank" else best["accession"]),
            sequence_source="refseq" if target["source"] == "refseq" else "genbank",
            gca_gcf_identity="verify_after_download" if gcf else ("ncbi_flags_different" if best["ncbi_flags_pair_different"] else "no_refseq_pair"),
            organism=target["organism"], tax_id=target["tax_id"], strain=target["strain"], bioproject=best["bioproject"],
            assembly_level=target["assembly_level"], records=target["records"], total_bp=target["total_bp"],
            annotation_provider=target["annotation_provider"], annotation_name=target["annotation_name"],
            atypical=best["atypical"] or target["atypical"], sequencing_tech=best["sequencing_tech"],
            linked_read_types=";".join(rtypes), n_usable_runs=sum(1 for l in links if l["biosample"] == bs and l["run_usable"]
                                                                  and l["assembly_accession"] == best["accession"]),
            other_assemblies_same_biosample=";".join(alts)))
        manifest.append(dict(accession=target["accession"], biosample=bs, purpose="primary",
                             include="genome,gff3,protein,seq-report"))
        if gcf:
            manifest.append(dict(accession=best["accession"], biosample=bs, purpose="gca_identity_check",
                                 include="genome,seq-report"))

    write_tsv(a.out / "assemblies.tsv.gz", sorted(asm.values(), key=lambda r: r["accession"]))
    write_tsv(a.out / "runs.tsv", run_rows)
    write_tsv(a.out / "assembly_run_links.tsv", links)
    write_tsv(a.out / "pool.tsv", pool_rows)
    write_tsv(a.out / "download_manifest.tsv", manifest)

    lk = [p for p in pool_rows if p["role"] == "linked"]
    summary = dict(
        schema="fna_pool_s0_v1", taxon=a.taxon, provenance=json.loads((raw / "provenance.json").read_text()),
        catalog=dict(refseq=sum(r["source"] == "refseq" for r in asm.values()),
                     genbank=sum(r["source"] == "genbank" for r in asm.values()),
                     genbank_biosamples=len(by_bs)),
        runs=dict(total=len(run_rows), by_platform=dict(Counter(r["instrument_platform"] for r in run_rows)),
                  usable=sum(r["usable"] for r in run_rows),
                  usable_by_read_type=dict(Counter(r["read_type"] for r in run_rows if r["usable"])),
                  usable_with_same_biosample_assembly=sum(1 for r in run_rows if r["usable"] and r["n_assemblies_same_biosample"]),
                  unusable_reasons=dict(Counter(x for r in run_rows for x in r["unusable_reasons"].split(";") if x))),
        links=dict(total=len(links), genus_disagrees=sum(not l["genus_agrees"] for l in links)),
        pool=dict(n_assemblies=len(pool_rows), n_unique_samples=len({p["biosample"] for p in pool_rows}),
                  by_role=dict(Counter(p["role"] for p in pool_rows)),
                  n_annotated_samples=sum(1 for p in pool_rows if p["annotation_provider"]),
                  n_linked_ONT_samples=sum("ont" in p["linked_read_types"].split(";") for p in lk),
                  n_linked_HiFi_samples=sum("pacbio_hifi_ccs" in p["linked_read_types"].split(";") for p in lk),
                  n_linked_CLR_samples=sum("pacbio_clr_subreads" in p["linked_read_types"].split(";") for p in lk),
                  n_linked_PacBio_unresolved_samples=sum("pacbio_unresolved" in p["linked_read_types"].split(";") for p in lk),
                  n_projects=len({p["bioproject"] for p in pool_rows}),
                  sequence_source=dict(Counter(p["sequence_source"] for p in pool_rows)),
                  assembly_level=dict(Counter(p["assembly_level"] for p in pool_rows)),
                  linked_assembly_level=dict(Counter(p["assembly_level"] for p in lk)),
                  atypical=sum(1 for p in pool_rows if p["atypical"])),
        download_manifest=dict(entries=len(manifest), by_purpose=dict(Counter(m["purpose"] for m in manifest))),
        caveats=["link_basis is exact BioSample only; whether the run was used for the assembly is not checked",
                 "GCA/GCF identity is checked by sequence after download, not inferred from NCBI pairing",
                 "no clinical/host screen, species exclusion or clone collapse at this stage"])
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "provenance"}, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["fetch", "build"])
    p.add_argument("--taxon", type=int, required=True)
    p.add_argument("--work", type=Path, required=True, help="cluster work dir; raw responses in <work>/raw")
    p.add_argument("--out", type=Path, help="repo dir for tables (build)")
    a = p.parse_args()
    fetch(a) if a.mode == "fetch" else build(a)


if __name__ == "__main__":
    main()
