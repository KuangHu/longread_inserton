# Next 10-sample call-out batch audit

Date: 2026-10-24.  Scope: freeze the immediate batch after the five-strain `nearby_competition_v2` acceptance run.  These files fix metadata only; raw FASTQ/FNA downloads still have to be verified on the cluster with the existing MD5/SHA256 checks.

## DD10

DD10 stays in the batch.  Its BioSample uses a MISAG/single-amplified template, but the 2026 genome announcement describes *Paraburkholderia* sp. strain DD10 as a cultured hyphospheric isolate from meadow-soil trapping of *Rhizophagus irregularis* mycelium, with Nanopore sequencing from LB liquid culture DNA.  The article links the same identifiers used here: SRA study SRP575236, BioSample SAMN47339132 and assembly GCF_049177125.1.

Action taken: `next_batch_proposal.tsv` preserves DD10 as an environmental-origin cultured isolate and records the metadata conflict as resolved by DOI `10.1128/mra.00286-25`.  We still keep the raw BioSample model in provenance; the resolution only decides whether DD10 can be run in the assembly+matched-ONT batch.

## Slot 10

The fixed rule was: choose the lowest-accession Burkholderiaceae genomic ONT WGS run with an exact BioSample-matched Complete Genome assembly, environmental or plant-associated cultured-isolate metadata, and a non-PCR library, after excluding the five acceptance strains and the nine already preselected batch strains.

The ENA query over taxon 119060 returned 867 genomic ONT runs from 817 BioSamples.  Joining exact BioSamples to NCBI GenBank Complete Genome assemblies gave 337 run/assembly rows from 316 BioSamples.  The first row that passed the rule was **Burkholderia seminalis 869T2**:

| field | value |
|---|---|
| sample_id | Bs869T2 |
| assembly | GCA_000705535.2 / GCF_000705535.2 |
| BioSample | SAMN02720836 |
| run | SRR14075918 |
| source | surface-sterilized vetiver roots |
| library | ONT MinION, WGS, GENOMIC, RANDOM |
| read bases | 522,326,261 |
| FASTQ MD5 | cf543453fcecc9329c8cc7c57e802c00 |

The complete genome has 3 records and a largest record of 3,012,498 bp.  It uses a RANDOM ONT library under the exact same BioSample as the assembly, so it satisfies the slot rule and is independent of the earlier *Paraburkholderia* pilot.

## Frozen files

| file | purpose |
|---|---|
| `next_batch_proposal.tsv` | the human-readable 10-row batch with caveats |
| `next_batch_inputs.tsv` | run metadata needed by download/verify jobs |
| `next_batch_selection_audit.tsv` | every exact BioSample match from the family-level Complete Genome scan, with the selection or exclusion reason |
| `next_batch_metadata_audit.json` | compact counts and the DD10 / slot-10 resolution |

This batch is source/data fixed.  The `input_status` column remains `metadata_frozen_download_verify_on_cluster` until `fetch_long_reads.py fetch` and the assembly downloader verify each file on the cluster.
