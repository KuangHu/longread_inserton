#!/usr/bin/env python3
"""Bounded real-read download probes. Prefix samples are NOT random samples.

ENA metadata resolves URLs; HTTPS ranges fetch a source prefix; complete read
records are rewritten to a new, complete gzip FASTQ. Source MD5 is deliberately
not claimed verified for prefix mode. Full mode downloads one complete FASTQ
with disk/size checks and verifies its archive MD5 before atomic publication.
"""
import argparse
import collections
import concurrent.futures
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import time
import urllib.parse
import urllib.request


SOURCES = [
    dict(id="marine_ccs", run="SRR28395891", platform="hifi", project="PRJNA1088973",
         evidence="https://doi.org/10.1186/s40793-025-00773-8"),
    dict(id="soil_ont", run="ERR15289757", platform="ont", project="PRJEB88618",
         evidence="https://doi.org/10.1038/s41467-026-69760-y"),
    dict(id="soil_hifi", run="ERR15289804", platform="hifi", project="PRJEB88618",
         evidence="https://doi.org/10.1038/s41467-026-69760-y"),
    dict(id="matched_gut_ont", run="ERR15285694", platform="ont", project="PRJEB88618",
         evidence="https://doi.org/10.1038/s41467-026-69760-y"),
    dict(id="matched_gut_hifi", run="ERR15289675", platform="hifi", project="PRJEB88618",
         evidence="https://doi.org/10.1038/s41467-026-69760-y"),
    dict(id="gut_102mags_hifi", run="SRR17687125", platform="hifi", project="PRJNA798244",
         evidence="https://doi.org/10.1038/s41467-022-34149-0"),
    dict(id="longitudinal_gut_ont", run="ERR14838512", platform="ont", project="PRJEB88320",
         evidence="https://doi.org/10.1038/s41586-025-09786-2"),
    dict(id="compost_ont", platform="ont", project="ONT_Compost_2026",
         url="https://ont-open-data.s3.amazonaws.com/compost_mgx_2026.04/basecalls/SUP/PAW97819.sup_v5.0.0.bam",
         source_size=126937137871, source_format="BAM",
         evidence="https://epi2me.nanoporetech.com/compost_mgx_2026.04/"),
]


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def request(url, headers=None):
    return urllib.request.urlopen(urllib.request.Request(url, headers=headers or {}), timeout=30)


def resolve(source):
    if "url" in source:
        return dict(source)
    fields = "run_accession,sample_accession,scientific_name,instrument_platform,instrument_model,library_source,library_strategy,base_count,read_count,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_format"
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?"+urllib.parse.urlencode(
        dict(accession=source["run"], result="read_run", fields=fields, format="json"))
    with request(url) as response:
        rows = json.load(response)
    rows = [r for r in rows if r["run_accession"] == source["run"]]
    if len(rows) != 1:
        raise ValueError("Run metadata missing or ambiguous")
    row = rows[0]
    paths = row["fastq_ftp"].split(";")
    sizes = row["fastq_bytes"].split(";")
    hashes = row["fastq_md5"].split(";")
    if len(paths) != 1 or not paths[0] or not (len(paths) == len(sizes) == len(hashes)):
        raise ValueError("This probe requires one archive FASTQ per selected run; inspect multi-file runs explicitly")
    return dict(source, url="https://"+paths[0], source_format="FASTQ.GZ",
                source_size=int(sizes[0]), archive_md5=hashes[0], ena_metadata=row,
                metadata_url=url, metadata_retrieved_utc=timestamp())


def fetch_range(url, start, end):
    """Require range semantics; never accidentally download a huge full file."""
    headers = {"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity"}
    with request(url, headers) as response:
        if response.status != 206:
            raise ValueError(f"Server did not honor Range: HTTP {response.status}")
        content_range = response.headers.get("Content-Range", "")
        if not content_range.startswith(f"bytes {start}-"):
            raise ValueError(f"Unexpected Content-Range: {content_range}")
        data = response.read(end-start+2)
        if len(data) != end-start+1:
            raise ValueError(f"Range size mismatch: {len(data)} != {end-start+1}")
        return data, {"status": response.status, "content_range": content_range,
                      "etag": response.headers.get("ETag"), "bytes": len(data)}


def read_fastq(path, limit=None):
    with gzip.open(path, "rt", encoding="ascii") as fh:
        n = 0
        while limit is None or n < limit:
            header = fh.readline()
            if not header:
                return
            seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
            if not qual:
                raise EOFError("Incomplete FASTQ record")
            header, seq, plus, qual = [s.rstrip("\r\n") for s in (header, seq, plus, qual)]
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual) or not seq:
                raise ValueError("Invalid four-line FASTQ record")
            if any(ord(c) < 33 or ord(c) > 126 for c in qual):
                raise ValueError("Quality characters outside Phred+33 range")
            yield header, seq, qual
            n += 1


def extract(prefix, kind, n):
    if kind == "FASTQ.GZ":
        records = list(read_fastq(prefix, n))
    else:
        import pysam
        records = []
        # The source is explicitly an incomplete BGZF prefix. A warning about
        # its absent terminal BGZF block is expected, never evidence of a full BAM.
        with pysam.AlignmentFile(str(prefix), "rb", check_sq=False, ignore_truncation=True) as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue
                seq, qualities = read.get_forward_sequence(), read.get_forward_qualities()
                if not seq or qualities is None:
                    continue
                records.append(("@"+read.query_name, seq, "".join(chr(q+33) for q in qualities)))
                if len(records) == n:
                    break
    if len(records) != n:
        raise EOFError(f"Prefix only contains {len(records)} complete records, need {n}")
    return records


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def qc(records):
    lengths = [len(s) for _, s, _ in records]
    qs = collections.Counter(ord(c)-33 for _, _, q in records for c in q)
    total = sum(lengths)
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= total/2:
            n50 = length
            break
    read_q = [-10*math.log10(sum(10**(-(ord(c)-33)/10) for c in q)/len(q)) for _, _, q in records]
    return dict(reads=len(records), bases=total, min_length=min(lengths),
                median_length=statistics.median(lengths), mean_length=total/len(records),
                max_length=max(lengths), subset_n50=n50,
                reads_ge_5kb=sum(n >= 5000 for n in lengths), reads_ge_10kb=sum(n >= 10000 for n in lengths),
                base_q_min=min(qs), base_q_max=max(qs), quality_values=len(qs),
                quality_histogram=dict(sorted(qs.items())),
                median_read_q_from_reported_qualities=statistics.median(read_q),
                all_qualities_identical=len(qs)==1,
                reads_with_ccs_in_header=sum("ccs" in h.lower() for h, _, _ in records),
                unique_read_ids=len({h.split()[0] for h, _, _ in records}),
                non_acgtn_bases=sum(c.upper() not in "ACGTN" for _, s, _ in records for c in s),
                first_read_headers=[h for h, _, _ in records[:3]])


def probe(source, args):
    target = args.out/source["id"]
    target.mkdir(parents=True, exist_ok=True)
    manifest = dict(source=source, status="running", started_utc=timestamp(),
                    mode="first_complete_records", requested_reads=args.reads,
                    sampling_warning="File-prefix sample, NOT random or representative of the run",
                    original_file_md5_verified=False, transfers=[])
    def save():
        (target/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    save()
    begin = time.monotonic()
    try:
        resolved = resolve(source)
        manifest["source"] = resolved
        save()
        prefix = target/"source.prefix"
        # Every invocation starts fresh; no stale prefix reused across remote versions.
        prefix.write_bytes(b"")
        maximum = min(resolved["source_size"], args.max_prefix_mb*1024*1024)
        records = None
        while prefix.stat().st_size < maximum:
            start = prefix.stat().st_size
            end = min(start+args.chunk_mb*1024*1024, maximum)-1
            for attempt in range(3):
                try:
                    data, transfer = fetch_range(resolved["url"], start, end)
                    break
                except Exception:
                    if attempt == 2:
                        raise
            if manifest["transfers"] and transfer["etag"] != manifest["transfers"][0]["etag"]:
                raise ValueError("Remote ETag changed between range requests")
            with prefix.open("ab") as fh:
                fh.write(data)
            manifest["transfers"].append(transfer)
            save()
            try:
                records = extract(prefix, resolved["source_format"], args.reads)
                break
            except (EOFError, OSError) as exc:
                manifest["last_incomplete_prefix_message"] = str(exc)
        if records is None:
            raise ValueError("Prefix byte cap reached before requested complete records")
        output = target/(source["id"]+".head"+str(args.reads)+".fastq.gz")
        temporary = output.with_suffix(output.suffix+".part")
        with gzip.open(temporary, "wt", encoding="ascii") as fh:
            for h, s, q in records:
                fh.write(f"{h}\n{s}\n+\n{q}\n")
        # Read to EOF: validates the newly created gzip CRC and all FASTQ records.
        checked = list(read_fastq(temporary))
        if checked != records:
            raise ValueError("FASTQ round-trip mismatch")
        temporary.replace(output)
        manifest.update(status="complete", completed_utc=timestamp(), qc=qc(checked),
                        downloaded_prefix_bytes=prefix.stat().st_size, output_file=output.name,
                        output_bytes=output.stat().st_size, output_sha256=digest(output),
                        output_gzip_crc_verified=True, seconds=time.monotonic()-begin)
        prefix.unlink()
        save()
    except Exception as exc:
        manifest.update(status="failed", error=f"{type(exc).__name__}: {exc}", seconds=time.monotonic()-begin)
        save()
    print(json.dumps({"id":source["id"],"status":manifest["status"],"error":manifest.get("error"),
                      "qc":manifest.get("qc"),"seconds":manifest["seconds"]}),flush=True)
    return manifest


def full_download(source, args):
    """One complete archive FASTQ, atomic and MD5 verified. No implicit BAM conversion."""
    resolved = resolve(source)
    if resolved["source_format"] != "FASTQ.GZ":
        raise ValueError("Full BAM conversion is not included; use the documented samtools workflow")
    size = resolved["source_size"]
    args.out.mkdir(parents=True, exist_ok=True)
    if size > args.max_full_gb*1e9 or shutil.disk_usage(args.out).free < size*1.1:
        raise ValueError("Full download exceeds configured size limit or available disk")
    output = args.out/Path(urllib.parse.urlsplit(resolved["url"]).path).name
    if output.exists():
        raise ValueError("Destination already exists")
    temporary = output.with_suffix(output.suffix+".part")
    h = hashlib.md5()
    received = 0
    with request(resolved["url"]) as response, temporary.open("wb") as fh:
        while block := response.read(1024*1024):
            received += len(block)
            if received > size:
                raise ValueError("Remote body exceeded metadata file size")
            fh.write(block)
            h.update(block)
    if received != size or h.hexdigest() != resolved["archive_md5"]:
        raise ValueError("Full file size or MD5 verification failed; file remains .part")
    output_with_metadata = dict(source=resolved, original_file_md5_verified=True,
                                bytes=received, retrieved_utc=timestamp())
    temporary.replace(output)
    output.with_suffix(output.suffix+".json").write_text(json.dumps(output_with_metadata,indent=2)+"\n")
    print(json.dumps(output_with_metadata),flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["probe","full"])
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--source",nargs="+",choices=[s["id"] for s in SOURCES])
    p.add_argument("--reads",type=int,default=500)
    p.add_argument("--workers",type=int,default=3)
    p.add_argument("--chunk-mb",type=int,default=4)
    p.add_argument("--max-prefix-mb",type=int,default=32)
    p.add_argument("--max-full-gb",type=float,default=10)
    args=p.parse_args()
    if min(args.reads,args.workers,args.chunk_mb,args.max_prefix_mb,args.max_full_gb)<=0:
        p.error("Limits must be positive")
    sources=[s for s in SOURCES if args.source is None or s["id"] in args.source]
    if args.mode=="full":
        if len(sources)!=1:p.error("Full mode requires exactly one --source")
        full_download(sources[0],args)
    else:
        if args.out.exists() and any(args.out.iterdir()):
            p.error("Probe output directory must be new or empty; stale results are not reused")
        args.out.mkdir(parents=True,exist_ok=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            results=list(pool.map(lambda s:probe(s,args),sources))
        (args.out/"summary.json").write_text(json.dumps(results,indent=2)+"\n")
        if any(r["status"]!="complete" for r in results):
            raise SystemExit(1)


if __name__=="__main__":
    main()

