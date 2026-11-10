"""Metadata-only inputs for the 50-sample panel: ENA ONT runs, exact-BioSample NCBI assemblies, FASTA stats.

Stage `query`: ENA read_run search over Burkholderiaceae (tax_tree 119060), all Oxford Nanopore runs;
NCBI Datasets GenBank genome reports for the family, kept when the BioSample equals that of a GENOMIC-source run.
Stage `download`: GCA FASTA (and paired GCF when asked) for the listed accessions; per-record
length and sequence SHA256 so continuity and GCA/GCF identity come from the files, not the reports.
Raw responses stay in the work directory on scratch; no reads and no discovery are run here.
"""

import argparse
import hashlib
import io
import json
import subprocess
import time
import urllib.parse
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
ENA_QUERY = 'tax_tree(119060) AND instrument_platform="OXFORD_NANOPORE"'
ENA_FIELDS = ("run_accession,study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,"
              "tax_id,scientific_name,strain,isolate,sample_title,instrument_platform,instrument_model,"
              "library_strategy,library_source,library_selection,library_layout,read_count,base_count,"
              "fastq_bytes,fastq_md5,fastq_ftp,submitted_bytes,host,host_scientific_name,isolation_source,"
              "environmental_medium,environment_biome,broad_scale_environmental_context,country,collection_date,first_public")


def curl(url, out=None):
    """GET with curl (urllib produced corrupted bodies on this cluster); retry three times."""
    for attempt in range(3):
        r = subprocess.run(["curl", "-sS", "--fail", "--max-time", "600", url] + (["-o", str(out)] if out else []),
                           capture_output=True)
        if r.returncode == 0:
            return r.stdout
        time.sleep(5 * (attempt + 1))
    raise SystemExit("curl failed: %s\n%s" % (url, r.stderr.decode()))


def query(work):
    raw = work / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    ena_url = ENA + "?" + urllib.parse.urlencode(dict(result="read_run", query=ENA_QUERY, fields=ENA_FIELDS,
                                                      format="tsv", limit=0))
    curl(ena_url, raw / "ena_ont_runs.tsv")
    lines = (raw / "ena_ont_runs.tsv").read_text().splitlines()
    head = lines[0].split("\t")
    runs = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    biosamples = sorted({r["sample_accession"] for r in runs if r["library_source"] == "GENOMIC" and r["sample_accession"]})
    # NCBI's BioSample genome endpoint returns nothing for SAMEA/SAMD accessions, so scan the whole
    # family (GenBank side; paired GCF is in each report) and join on the exact BioSample accession.
    wanted, reports, pages, token = set(biosamples), [], [], None
    while True:
        url = NCBI + "/taxon/119060/dataset_report?" + urllib.parse.urlencode(
            dict(page_size=1000, **{"filters.assembly_source": "genbank"}, **({"page_token": token} if token else {})))
        d = json.loads(curl(url))
        page = d.get("reports", [])
        pages.append(dict(page=len(pages) + 1, reports=len(page), total_count=d.get("total_count")))
        reports += [r for r in page if r["assembly_info"].get("biosample", {}).get("accession") in wanted]
        token = d.get("next_page_token")
        if not token:
            break
    scanned = sum(q["reports"] for q in pages)
    if scanned != pages[0]["total_count"]:
        raise SystemExit("NCBI scan incomplete: %d of %s" % (scanned, pages[0]["total_count"]))
    with open(raw / "ncbi_biosample_reports.jsonl", "w") as h:
        for r in reports:
            h.write(json.dumps(r) + "\n")
    prov = dict(date=time.strftime("%Y-%m-%d"), ena_url=ena_url, ena_query=ENA_QUERY, ena_rows=len(runs),
                ena_genomic_biosamples=len(biosamples),
                ncbi_scan=NCBI + "/taxon/119060/dataset_report?page_size=1000&filters.assembly_source=genbank (all pages)",
                ncbi_pages=pages, ncbi_genbank_reports_scanned=scanned, ncbi_reports_exact_biosample=len(reports),
                ena_tsv_sha256=hashlib.sha256((raw / "ena_ont_runs.tsv").read_bytes()).hexdigest())
    (raw / "query_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps({k: v for k, v in prov.items() if k != "ncbi_pages"}, indent=2))


def fasta_stats(data):
    recs, name, seq = [], None, []
    for line in data.decode().splitlines():
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        recs.append((name, "".join(seq)))
    lens = sorted((len(s) for _, s in recs), reverse=True)
    total, acc, n50 = sum(lens), 0, 0
    for L in lens:
        acc += L
        if acc * 2 >= total:
            n50 = L
            break
    return dict(records=len(recs), total_bp=total, largest_bp=lens[0] if lens else 0, n50_bp=n50,
                n_bases=sum(s.upper().count("N") for _, s in recs),
                record_sha256=sorted(hashlib.sha256(s.upper().encode()).hexdigest() for _, s in recs),
                record_lengths={n: len(s) for n, s in recs})


def download_one(acc, out):
    target = out / (acc + ".fna")
    meta = out / (acc + ".json")
    if meta.exists():
        return json.loads(meta.read_text())
    url = NCBI + "/accession/" + acc + "/download?include_annotation_type=GENOME_FASTA"
    data = curl(url)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        bad = z.testzip()
        files = [n for n in z.namelist() if n.endswith("_genomic.fna") and acc in n]
        if bad or len(files) != 1:
            return dict(accession=acc, status="download_failed", error="crc %s members %s" % (bad, files))
        seq = z.read(files[0])
    target.write_bytes(seq)
    res = dict(accession=acc, status="downloaded", url=url, fasta=str(target),
               fasta_sha256=hashlib.sha256(seq).hexdigest(), **fasta_stats(seq))
    meta.write_text(json.dumps(res) + "\n")
    return res


def download(work, accessions):
    out = work / "genomes"
    out.mkdir(parents=True, exist_ok=True)
    accs = [a for a in Path(accessions).read_text().split() if a]
    with ThreadPoolExecutor(max_workers=4) as pool:
        res = list(pool.map(lambda a: download_one(a, out), accs))
    print(sum(r["status"] == "downloaded" for r in res), "downloaded of", len(accs))
    for r in res:
        if r["status"] != "downloaded":
            print(r)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["query", "download"])
    p.add_argument("--work", type=Path, default=Path("/global/scratch/users/kh36969/longread_callout/batch50_meta"))
    p.add_argument("--accessions", help="download: file of GCA/GCF accessions")
    a = p.parse_args()
    query(a.work) if a.stage == "query" else download(a.work, a.accessions)


if __name__ == "__main__":
    main()
