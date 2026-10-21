# Real CCS pilot: SRR28395891 vs marine study MAG pool

This is an exploratory real-data pipeline run, **not a fully verified same-sample assembly/MAG benchmark**. The source is a natural Mediterranean marine microbial community, not a laboratory strain mixture.

## Provenance limitation

The intended design was same-sample assembly/MAG plus HiFi/CCS reads. Metadata inspection found that the MAG BioSample `derived-from` field can list all four run accessions together, and the paper's MAG table does not provide a unique source-sample column. We therefore cannot establish that each reference MAG was assembled specifically from SRR28395891. The pipeline uses the paper's **study-wide MAG pool** as a competitive reference and must not be described as strict same-sample validation. References can come from other samples in the study, and this affects interpretation of allele differences.

The alternative 2019 project PRJNA674982 was checked but returned no assemblies in the ENA and NCBI assembly queries used here; that does not prove the sequences do not exist elsewhere. No 2019 read run was used in this pilot.

## Inputs

- Reads: SRR28395891, Med-OCT2021-15m, PRJNA1088973. Whole archived FASTQ, not a prefix subsample. Paper describes CCS processing; `_subreads.fastq.gz` in the archive filename is not used to infer CLR.
- Reference: Zenodo record 14650439, `MAG-LRMED.tar`, 605,245,952 bytes, MD5 `51d2fd500477514ce776374983b974c1`.
- The archive contains 419 FASTA files. Only the 418 MAGs present in paper Table S2 are used. `Mycoplasmatales-LRMED-162` is excluded because it is absent from that table.
- Reference: 17,441 contigs, 594,150,817 bp. Contig names are prefixed with MAG names to avoid collisions such as repeated `contig_1`.
- Reference SHA256: `adf45f5c5e21d5e30b847bf3fb7ea6cffe63c90907e0e9ca9b661a49cf1b448b`.
- Archive vs Table S2 validation found 14 MAGs with a length or contig-count discrepancy. The archived sequences remain in the competitive reference; all reference lengths and coverage denominators come from the actual FASTA, not the table. See `table_validation.json`; events on these MAGs need an additional provenance flag.
- This is a MAG collection, not the complete original metagenome assembly; unbinned sequences are absent.

Sources:
- https://doi.org/10.1186/s40793-025-00773-8
- https://zenodo.org/records/14650439
- https://mgegroup.com.es/datasets/
- https://www.ebi.ac.uk/ena/browser/view/SRR28395891

## Reproduction

From the repository root, with the pipeline dependencies installed:

```bash
python validation/real_hifi_pilot/fetch_reference.py /path/to/hifi_marine_pilot
python validation/real_hifi_pilot/prepare_reference.py /path/to/hifi_marine_pilot
python tools/downloads/download_probe.py full --source marine_ccs \
  --out /path/to/hifi_marine_pilot/reads --max-full-gb 1
python validation/real_hifi_pilot/align.py /path/to/hifi_marine_pilot /path/to/minimap2
python validation/real_hifi_pilot/summarize.py /path/to/hifi_marine_pilot --coverage-only
python validation/real_hifi_pilot/call_targets.py /path/to/hifi_marine_pilot \
  /path/to/sniffles /path/to/insertion-site-recovery
python validation/real_hifi_pilot/summarize.py /path/to/hifi_marine_pilot
```

The initial sequential FASTQ transfer timed out after 413,138,944 bytes. `resume_reads.py` records the bounded HTTP Range recovery used here: three parallel workers, 8 MiB ranges, exact Content-Range and byte-count checks, followed by MD5 of the entire assembled file before publication. A partial file is never passed to the caller.

`fetch_reference.py` additionally requires openpyxl. `run.sh` reproduces the final staged workflow using this workspace layout. `run_initial_attempt.sh` records the monolithic attempt that was stopped after alignment and during broad calling. The successful targeted calling and recovery commands are recorded in results/calling_manifest.json. Caller thresholds are unchanged from prototype defaults. The initial minimap2 run was interrupted as memory approached the 8 GiB environment limit; no OOM kill was recorded. The completed alignment used the explicit wrapper `minimap2_lowmem.sh` to set `-K 50M` (50 million query bases per batch), keeping the full reference in one index. This changes read batching rather than alignment thresholds. A fresh output directory is required. Downloads, FASTQ, reference and BAM are not committed to Git.

Coverage summaries use primary, nonduplicate, non-QC-failed alignments with MAPQ >=20. Aligned M/=/X blocks count toward depth; deletions and skipped regions do not. Mean depth across a MAG is not a guarantee of coverage at a candidate junction. The reported quality values in the CCS archive can be uniform; they are not an independent measurement of accuracy.

A `read_supported_pair` is computational local evidence for both structures. It is not proof of a mobile element, evolutionary direction, or same-strain/cell origin. Empty-site sequence is reference-derived, not a strain-resolved read consensus.

## Archive-specific observations

The complete FASTQ contains 418,396 reads and 1,764,906,237 bases (median 3,410 bp; N50 5,309 bp). Every header contains the authors’ `CCS15` label, and all reported quality characters are Q30. The paper Table S1 reports a different analyzed subset (121,917 reads for this sample); this pilot uses the full archived run and does not claim to reproduce the paper’s filtering or counts.

## Results

See [the Chinese report](../../docs/REAL_HIFI_PILOT.md) for results and limitations. (From this directory, the report is at `../../docs/REAL_HIFI_PILOT.md`.)

The target selection was frozen by coverage and assembly QC before targeted SV calling: top 10 eligible MAGs, 235 contigs. Final results: 312 SV records (276 INS/DEL); 20 sequence representations recovered; 3 strict read-supported pairs and 17 repeat-boundary ambiguous candidates, of which 13 also have sufficient support for both alleles.
