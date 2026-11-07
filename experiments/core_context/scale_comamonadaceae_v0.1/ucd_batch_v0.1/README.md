# U–X–D comparison of the existing candidate catalogue, with source and independent checks (ucd_batch v0.1)

**Scope:** every `match_both_lt2500` pair of the existing candidate catalogue (`match_both_lt2500_pairs_b1-2-3.tsv.gz`), 95,328 pairs.
- This is a full comparison of that catalogue. It is **not** a complete discovery over every FASTA.
- Pairs are taken as listed. Nothing is filtered by old read-check, gap-caller or S6 results; `s6_grouping` is a carried label.
- The catalogue itself comes from earlier steps (core match, boundary review, span class). Those steps were not rerun.

**Acceptance conditions:**
1. The extracted sequences correspond exactly to the source assembly.
2. The insertion is similar, and upstream and downstream both differ.

`final_status = pass` requires both. Results are a structural pattern: no empty site, no mobility, no read evidence.

## Run

```bash
D=scale_comamonadaceae_v0.1/ucd_batch_v0.1; W=/global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.1
sbatch --cpus-per-task=4 --mem=32G ucd_batch.sbatch freeze      # input table, source verification, chunks, hashes
sbatch --array=0-47 ucd_batch.sbatch chunk                       # resubmitting skips finished chunks
sbatch --cpus-per-task=4 --mem=48G ucd_batch.sbatch merge        # per-ID completeness, then tables
python ucd_batch.py cleanup --apply --out $D --work $W           # only after a complete, hash-checked merge
python ucd_batch_examples.py --out $D --work $W
```

**Frozen before the run (`freeze.json`):**
- the sha256 of `ucd_compare.py`, `ucd_checks.py` and `ucd_batch.py` (a chunk refuses to run if the code differs);
- CONFIG and CHECK_CONFIG;
- the catalogue and manifest sha256;
- the input pair table (`input_pairs.tsv.gz`) and the source table (`sources.tsv.gz`).

**Source files:** all 1,669 assembly FASTAs were checked against the manifest's verified sha256, and all are `verified`. Nothing was downloaded.

## What each row records

The tables have one row per pair and interval kind. The two kinds are kept apart in separate files:
- `cds_core` is the catalog CDS.
- `candidate_insertion` is the CDS plus the pair's own resolved extensions. It is pair-specific.

| block | columns |
|---|---|
| provenance per side (a, b) | occurrence, sample, assembly accession.version, FASTA sha256, contig, contig length, strand, U / X / D intervals, contig-edge flag, sha256 (16 hex) of the oriented U / X / D |
| screen (v0.1, kept as is) | `screen_verdict`, insertion identity and coverages, U/U and D/D, reasons |
| Check 1 | `source_check` (pass / fail / not_run), `source_reasons` |
| Check 2 | `ucd_check`, `chk_*`: insertion chain identity and coverage, union coverage, orientation, usable bp, the 8 context comparisons (class, identity, coverage of each flank), reasons |
| result | `final_status`, `failure_reason`, `screen_check_agree` |

**Coordinates** are 0-based, half-open, on the forward strand of the contig as in the source FASTA. U, X and D are given in the insertion orientation; a minus strand is the reverse complement of the forward interval.

## Checks

**Check 1, source and extraction** (`ucd_checks.Source`, `check1_side`):
- The FASTA sha256 is recomputed on every load and must equal the manifest.
- The file is indexed from its bytes. Every interval is read again by byte offset, and the minus strand is complemented with a full IUPAC table.
- The expected U / X / D intervals are recomputed from the input coordinates. Intervals and bases must equal the screen's, base by base.
- Errors are reported, never repaired. Recorded errors: checksum mismatch, missing contig, out-of-bounds or empty interval, irregular line widths, interval or sequence different from the screen.
- A contig edge that shortens a flank is recorded, not an error. The usable-context rule then judges that flank.

**Check 2, independent comparison** (`ucd_checks.check_ucd`): a separate implementation, run on the Check-1 sequences.
- **Insertion:** coverage must come from one collinear chain of reliable HSPs (≥ 90% identity, ≥ 30 bp) on one strand.
  - Each next HSP starts and ends later in both sequences. An overlap of at most 20 bp (`max_chain_overlap`) is allowed and counted once.
  - Repeated hits, out-of-order or inverted pieces are not added.
  - `insertion_similar` needs the chain to cover ≥ 90% of X and of X′.
  - If only the union of all pieces reaches 90%, the row is `rearranged_or_fragmented_insertion`.
- **Contexts:** the frozen v0.1 rule: 1 kb windows; union of HSPs ≥ 80% on each flank.
  - `different` means ≤ 20% on both flanks; `similar` means ≥ 80% on both; anything between is `uncertain`.
  - Comparisons: same side, inverted, and cross side on both strands.
  - Also as in v0.1: an insufficient usable flank (< 800 bp) and a shared physical flank.
- **Not required:** an empty site, mobility, or any function name.

## Integrity and timing

- **Completeness** (`completeness.json`): 190,656 of 190,656 pair_id × interval kind, each exactly once. No duplicates, nothing unexpected.
- **Runs** (`jobs.txt`, `chunks_done.json`):
  - freeze: 40 min (job 27603383);
  - chunk 0 smoke test: 5 min (27605962);
  - chunks 1–47: array 27606582;
  - merge: 10 min (27608682).
- **Failures and retries:** array task 42 failed on n0004 before doing any work: the scratch conda environment was not found on that node (`EnvironmentLocationNotFound`). It was retried as 27608032_42 with `--exclude=n0004.lr6` and completed on n0003. No other task failed. Methods, parameters and inputs were not changed.
- **Wall time:** first chunk start 12:14, last chunk done 12:53, merge done 13:06.
- **Compute:**
  - 4.5 node-hours of chunk tasks;
  - 172 allocated core-hours (whole lr6 nodes of 32–40 cores);
  - 28 CPU-hours actually used by the workers, mostly blastn. The rest is I/O and process start-up on scratch.
- **Memory:** ≤ 0.4 GB RSS per chunk task.
- **Cleanup** (`cleanup.json`): the chunk directories (98.6 MB, rebuildable) were removed after the merged files matched their sha256 in `summary.json`.

## Results (`summary.json`)

Counts are pairs, then distinct occurrences / samples / cores. Pairs are not independent events: one occurrence is in many pairs, and a two-gene core can repeat the locus of a single-gene core.

**Final table**, final_status = pass:

| | cds_core | candidate_insertion |
|---|---|---|
| **pass** | **48,411** / 20,572 / 1,589 / 3,383 | **56,811** / 23,829 / 1,616 / 4,099 |
| same assembly, single gene | 12,395 / 7,723 / 1,052 / 895 | 15,394 / 9,317 / 1,064 / 1,162 |
| same assembly, two gene | 892 / 685 / 227 / 150 | 1,030 / 845 / 251 / 201 |
| cross assembly, single gene | 31,097 / 11,556 / 1,522 / 2,229 | 35,579 / 12,864 / 1,553 / 2,543 |
| cross assembly, two gene | 4,027 / 2,775 / 1,051 / 501 | 4,808 / 3,450 / 1,139 / 672 |

**Every row** (`ucd_check`; all rows kept in scratch, see Files):

| ucd_check | cds_core pairs | candidate_insertion pairs |
|---|---:|---:|
| ucd_supported | 48,411 | 56,811 |
| insertion_dissimilar | 32,389 | 36,550 |
| uncertain_context | 11,200 | 1,096 |
| only_one_side_different | 2,572 | 43 |
| insufficient_context | 483 | 459 |
| insertion_partial | 73 | 160 |
| swapped_or_rearranged_context | 60 | 67 |
| same_context | 27 | 0 |
| rearranged_or_fragmented_insertion | 0 | 0 |
| not run: source_check failed | 113 | 142 |

**Source check failures** (113 cds_core / 142 candidate_insertion pairs):
- **Screen extraction differs** (52 occurrence sides with 117 U / X / D parts in cds_core; 54 sides with 121 parts in candidate_insertion; 26 assemblies).
  - Every side is on the minus strand.
  - Every differing part contains IUPAC codes (R, Y, K, M, S, W), as checked again from the source.
  - The v0.1 screen reverses those letters without complementing them; Check 1 complements them.
  - This is a defect of the screen's extraction, which Check 1 caught. These pairs are not passed.
- **interval_out_of_bounds** (25 sides, candidate_insertion only): the pair's extension reaches past a contig end, at the start or end of a contig, or on a short contig. Not repaired.

**Screen vs Check 2:**
- `cds_core`: identical on every row where both ran.
- `candidate_insertion`: 26 rows are `insertion_partial` in the screen and `insertion_dissimilar` in Check 2.
  - All 26 have a single ~30 bp HSP at exactly 90% identity covering 2.5% of X.
  - Check 2 searches A's three sequences against B's three in one call. The larger search space puts that HSP over the E-value cutoff.
  - None of the 26 passes either way.

**Why the two interval kinds differ:**
- On `cds_core`, 11,200 pairs are uncertain and 2,572 one-side-different.
- On `candidate_insertion`, these are 1,096 and 43: the shared sequence often continues past the CDS, so CDS flanks start inside it.
- A CDS is never called a complete insertion. Both kinds are reported side by side.

**insertion_dissimilar:** every `insertion_dissimilar` pair has core match `homolog_full_coverage`.
- That is 32,389 of the 78,998 such pairs on `cds_core`.
- These core copies are related at protein level but have no collinear ≥ 90% DNA match over ≥ 90%, so they fail condition 2.
- Of the 16,330 `identical_core_dna` pairs, none is dissimilar.
- These pairs are kept in the all-rows tables.

**Length note:** 251 of 190,656 candidate intervals (0.13%) are 2,500–2,510 bp. The catalogue's `lt2500` class is from the occurrence-level span estimate, while the interval here is the pair's own extension; they differ by a few bp. These rows were not removed (no refiltering of inputs).

## Examples (`examples.md`)

There is one `candidate_insertion` pass per pair scope × structure, chosen by a fixed rule: the longest insertion, then the smallest pair_id.

Each example shows:
- the provenance of both sides;
- identity and coverage of the insertion and all 8 context comparisons;
- the four junctions (40 bp each side);
- the Check-2 HSPs.

The rule picks long elements, so two examples (cross-assembly single and two gene) show the same 2,501 bp locus seen through two cores.

## Files

**In this directory (git):**
- `freeze.json`, `input_pairs.tsv.gz`, `sources.tsv.gz`
- `ucd_final_cds_core.tsv.gz`, `ucd_final_candidate_insertion.tsv.gz` (pass rows, all columns)
- `summary.json`, `completeness.json`, `chunks_done.json`, `jobs.txt`, `cleanup.json`, `examples.md`

**In scratch** (`/global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.1`; paths, bytes and sha256 in `summary.json`):

| file | content | size |
|---|---|---|
| `ucd_all_rows_cds_core.tsv.gz`, `ucd_all_rows_candidate_insertion.tsv.gz` | every row, including uncertain, failed and not run | 28 MB, 27 MB |
| `evidence_sequences.fa.gz` | Check-1 U / X / D of every compared occurrence interval, deduplicated; headers carry assembly, contig, interval, strand, FASTA sha256 and sequence sha256 | 52 MB |
| `evidence_hsps.tsv.gz` | every Check-2 HSP, with `in_insertion_chain` | 4.8 MB |

## Cost of the remaining full run

The boundary review has 3,894,917 core pairs, so **3,799,589 remain**.

| measure | per pair | remaining pairs |
|---|---|---|
| node time (32–40 cores) | 0.171 node-s | ≈ 180 node-h |
| allocated core time | 6.5 core-s | ≈ 6,800 core-h |
| CPU used | 1.1 CPU-s | ≈ 1,150 CPU-h |
| evidence + all-row tables | ≈ 1.2 KB | ≈ 4.5 GB |

- **Wall time:** about 4 h at 48 concurrent array tasks (≈ 800 chunks of 5,000 pairs).
- **This is an upper bound:**
  - most remaining pairs are `match_boundary_unknown_pending`, with no `candidate_insertion` interval, so only `cds_core` is compared for them;
  - batching more pairs per blastn call would reduce the start-up and I/O overhead (allocated core time is about 6× the CPU used).
- `freeze` would need the catalogue of that set. The rest of the pipeline is unchanged.
