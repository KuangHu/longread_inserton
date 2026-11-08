# U–X–D comparison of the whole boundary-review catalogue (ucd_batch v0.2)

**Supersedes:** `../ucd_batch_v0.1/` (commit faa63c3). That directory and its scratch results are kept unchanged as
history. For the 95,328 pairs of v0.1, this version reuses the v0.1 rows, except the 117 pairs whose sequences the
v0.1 screen extracted wrongly; those were rerun with the fixed code.

**Scope:** every core pair of the boundary review (b1-2-3), 3,894,917 pairs.
- This is a complete comparison of that existing catalogue. It does **not** claim to find every possible insertion
  in all FASTA.
- Nothing is filtered by old read-check, gap-caller or S6 results. Nothing infers a new boundary.

**Two acceptance conditions, and no others:**
1. The U / X / D sequences equal the stated intervals of the source assembly (Check 1).
2. The insertion is similar, and upstream and downstream both differ (Check 2).

There is no empty-site, read-support, gene-completeness, S6 or mobility condition.

## Order of work

| step | directory | what |
|---|---|---|
| 1. fix | `../../ucd_compare.py` v0.1.1 | `rc()` complements every IUPAC code (R↔Y, K↔M, B↔V, D↔H; S, W, N unchanged) and keeps case. v0.1 reversed ambiguity codes without complementing them. Coordinates unchanged. Regression tests: minus-strand extraction with ambiguity codes, rc(rc(s)) = s, A/C/G/T results unchanged (`test_ucd_compare.py`). Check 1 is unchanged. |
| 2. rerun | `rerun/` | the 117 pairs that failed Check 1 with `screen_sequence_differs` in faa63c3 (113 in both kinds, 4 in candidate_insertion only), both interval kinds, plus 2,000 random unaffected control pairs (seed 20261107) |
| 3. strict export of the first batch | `first_batch/` | v0.1 rows reused + the 117 rerun pairs; strict < 2,500 bp table |
| 4. remaining pairs | `remaining/` | the frozen catalogue minus the 95,328 processed IDs |
| 5. final | `final/` | all of the above combined, checked per pair_id × interval kind |

## 1–2. Fix and rerun (`rerun/rerun_report.json`, `rerun/rerun_affected.tsv`)

| | cds_core | candidate_insertion |
|---|---|---|
| affected pairs | 117 | 117 |
| source_check before → after | 113 fail, 4 pass → 117 pass | 117 fail → 117 pass |
| ucd_check after | 55 ucd_supported, 48 insertion_dissimilar (+4 already so), 8 uncertain_context, 2 only_one_side_different | 44 ucd_supported, 72 insertion_dissimilar, 1 uncertain_context |
| **new passes** | **55** | **44** |
| control pairs identical to faa63c3 in every compared column | 2,000 / 2,000 | 2,000 / 2,000 |

- Out-of-bounds pairs (25 in faa63c3) were not affected and stay failed. No coordinate was cut, moved or wrapped.
- Unaffected v0.1 rows are reused with `code_version = ucd_batch_v0.1@faa63c3`.
- New rows carry `ucd_batch_v0.2+48a91ae6ca23`: v0.2 plus the hash of the three frozen code files, sha256 in each
  `freeze.json`.

## 3. Strict length rule

A pair is in the strict table if:
- `source_check = pass`;
- `ucd_check = ucd_supported`;
- the X actually extracted on **both** sides has 0 < length < 2,500 bp.

Length is computed from the Check-1 interval and must equal the extracted sequence's length (all rows agree). The
occurrence-level `lt2500` label is not used, and exactly 2,500 bp fails.

A CDS / core under 2,500 bp does not mean the complete insertion is under 2,500 bp. The two interval kinds are saved
and counted separately.

**First batch** (95,328 pairs; `first_batch/combined_summary.json`):

| | cds_core | candidate_insertion |
|---|---|---|
| pass (pairs / occurrences / samples / cores) | 48,466 / 20,611 / 1,590 / 3,388 | 56,855 / 23,865 / 1,616 / 4,105 |
| strict < 2,500 | 48,466 / 20,611 / 1,590 / 3,388 | 56,719 / 23,865 / 1,616 / 4,105 |
| excluded: pairs / single-side intervals (distinct) | 0 / 0 | 136 / 250 (234) |

## 4. Remaining pairs (`remaining/freeze.json`, `remaining/merge.json`)

- **ID difference:** 3,894,917 catalogue pairs minus 95,328 processed IDs (all present in the catalogue) =
  **3,799,589**. There were no ID collisions and no input errors.
- **Sources:** all 1,678 assembly FASTAs were verified against the manifest sha256.
- **Pairs with a candidate interval:** both interval kinds were run.
- **Pairs without one:** only `cds_core` was compared. The `candidate_insertion` row is `not_run` and names the
  reason, e.g. `candidate_interval_not_available(boundary_up=inconsistent_between_partners,...)` or
  `pair_extension_not_determined`.

## 5. Final result (`final/combined_summary.json`, `final/verification.json`)

**Completeness** per pair_id × interval kind, for each kind:
- 3,894,917 expected, 3,894,917 observed;
- 0 missing, 0 duplicated, 0 unexpected;
- 0 rows without a status, and 0 rows whose length disagrees with the interval.

Counts are pairs / distinct occurrences / samples / cores. Pairs are not independent events: an occurrence sits in many pairs, and a two-gene core can repeat a single-gene locus.

**All records by ucd_check:**

| ucd_check | cds_core pairs | candidate_insertion pairs |
|---|---:|---:|
| ucd_supported (final pass) | 1,045,058 | 86,528 |
| insertion_dissimilar | 2,498,385 | 40,873 |
| uncertain_context | 220,880 | 1,250 |
| only_one_side_different | 80,945 | 54 |
| insertion_partial | 44,724 | 263 |
| insufficient_context | 2,858 | 473 |
| same_context | 1,507 | 0 |
| swapped_or_rearranged_context | 458 | 88 |
| rearranged_or_fragmented_insertion | 102 | 0 |
| not run: source_check failed | 0 | 54 (interval_out_of_bounds) |
| not run: candidate interval not available | — | 3,765,334 |
| **total** | **3,894,917** | **3,894,917** |

**Passing both checks, and the strict table:**

| | cds_core | candidate_insertion |
|---|---|---|
| **pass** | 1,045,058 / 110,789 / 1,669 / 7,805 | 86,528 / 28,100 / 1,621 / 4,714 |
| pass, same assembly | 28,062 / 14,068 / 1,275 / 1,560 | 23,602 / 12,511 / 1,088 / 1,658 |
| pass, cross assembly | 1,016,996 / 108,776 / 1,664 / 7,327 | 62,926 / 19,209 / 1,589 / 3,638 |
| **strict < 2,500 bp** | **1,045,058 / 110,789 / 1,669 / 7,805** | **57,145 / 24,194 / 1,617 / 4,177** |
| strict, same assembly, single gene | 24,168 / 11,956 / 1,222 / 1,307 | 15,388 / 9,390 / 1,068 / 1,178 |
| strict, same assembly, two gene | 3,894 / 2,112 / 599 / 253 | 1,006 / 859 / 254 / 201 |
| strict, cross assembly, single gene | 929,274 / 92,328 / 1,660 / 5,896 | 35,935 / 13,146 / 1,556 / 2,606 |
| strict, cross assembly, two gene | 87,722 / 16,448 / 1,555 / 1,431 | 4,816 / 3,467 / 1,144 / 677 |
| excluded by length: pairs / single-side intervals (distinct) | 0 / 0 | 29,383 / 58,721 (15,570) |

**Reading the two kinds:**
- `cds_core` pass means the coding core and its 1 kb flanks meet the rule. It is not a complete insertion, and most
  of these pairs have no resolved element boundary.
- The `candidate_insertion` strict table is the target set: a resolved element under 2,500 bp on both sides,
  similar, with different flanks.
- The remaining run had 34,255 pairs with a candidate interval. 29,673 of them pass, and 29,247 of those are
  excluded by length: at least one side's element is ≥ 2,500 bp.

**Unresolved categories** (kept with reasons in the all-rows tables; none is counted as different context):
- `uncertain_context` (between thresholds, or a cross / inverted comparison uncertain)
- `insufficient_context`
- `insertion_partial`
- `rearranged_or_fragmented_insertion`
- `interval_out_of_bounds` (54 candidate pairs)
- `candidate_interval_not_available` (3,765,334 pairs on the candidate kind; their `cds_core` rows are complete)

**Screen vs Check 2** (`final/verification.json`):
- `cds_core`, 104 rows:
  - 90 are `ucd_supported` in the screen but `rearranged_or_fragmented_insertion` in Check 2. The screen's union
    coverage counts out-of-order pieces; Check 2's collinear chain does not, so these pairs do not pass.
  - 9 + 3 more rows are rearranged in Check 2.
  - 2 are uncertain → ucd_supported. Check 2's combined search drops borderline cross HSPs.
- `candidate_insertion`: 26 rows, partial → dissimilar (as in v0.1).
- Final status follows Check 2.

## Examples (`final/examples.md`)

Four strict `candidate_insertion` pairs, one per scope × structure, chosen by fixed rule (the longest insertion, then
the smallest pair_id). Each shows:
- assembly accession.version, contig, coordinates, strand and FASTA sha256 of both sides;
- X/X′ identity and coverage (one collinear HSP, 95.3–100%, 2,498–2,499 bp) and all 8 context comparisons;
- the four junctions;
- the Check-2 HSPs.

## Long-read evidence of the strict pairs (`final/read_evidence_summary.json`, `final/strict_candidate_read_evidence.tsv.gz`)

`../../ucd_read_evidence.py` adds evidence columns to the 57,145 strict `candidate_insertion` pairs. Pass / fail is
unchanged (57,145 / 57,145 still `pass`). Long-read support is **not** an admission condition. It is reported apart
from the assembly result.

**Interval-level read check, round 5 v0.2** (`rc_<side>_*`):
- **Matching rule:** a side takes a result only when the tested interval is this row's X: the same sample (same
  manifest assembly), contig, lo / hi and occurrence, with the strand checked against round 5. No strand mismatch
  occurred.
- **No inheritance:** if the occurrence was tested on a different interval, nothing is inherited. Those intervals
  are listed in `rc_<side>_other_intervals`.
- **Copied as they are:** the raw interval class, verdict, count mode (molecules / read records / assumed
  molecules), P_matched, A_matched and locus_ambiguous.

| side status | sides |
|---|---:|
| supported | 29,018 |
| insufficient_coverage | 269 |
| locus_ambiguous | 34 |
| tested_other_interval | 46 |
| not_tested (the assembly has usable long-read runs) | 119 |
| no_reads (no usable long-read run for the BioSample) | 84,804 |

- **Pairs:**
  - 7,803 with both sides supported;
  - 13,264 with one side supported and the other `no_reads`;
  - 35,655 with `no_reads` on both sides.
  - The full breakdown, also per scope, is in the summary.
- No conflict occurred among the strict pairs.

**Read-site search, v0.3.2** (`rss_<side>_*`), as a separate field:
- **Coverage:** the 133 real samples (validation 1–2, scale 1–3); pilot samples are labelled development.
- **Rule:** flank 1 kb, at least 3 molecules, a `known_catalog_position*` position listing this occurrence.

| side status | sides |
|---|---:|
| recovered | 12,892 |
| not_recovered (not counted as unsupported) | 2,662 |
| seen_below_threshold | 97 |
| sample_not_searched | 98,639 |

## Integrity, retries and cost (`jobs.txt`, `*/chunks_done.json`, `*/merge.json`)

- **Jobs:**
  - freeze, rerun set: 4.6 min (27615001);
  - freeze, remaining set: 11 min, 7.6 GB RSS (27615002);
  - rerun: 3 chunks, 1.5–2.5 min each (27615261);
  - remaining: 760 chunks of 5,000 pairs (27615539, 64 at a time, 15:54–18:36);
  - merge: 31 min (27622658);
  - combine: 13 min (27622659).
- **Failures:** none in v0.2. All 767 jobs and tasks completed on the first attempt, so nothing was retried. The v0.1
  retry (task 42, node n0004) is recorded in `../ucd_batch_v0.1/`.
- **Measured cost of the remaining run** (3,799,589 pairs):
  - 147.5 node-hours, 5,755 allocated core-hours, 745 CPU-hours used by the workers;
  - 0.14 node-s per pair.
  - Most of these pairs ran `cds_core` only.
- **Frozen code:** `ucd_compare.py` 4909ad16…, `ucd_checks.py` 35111726…, `ucd_batch.py` a58aec7f…. Same for the
  rerun and remaining runs; a chunk refuses to run on other code. Committed in 9ef3efb before the remaining run.

## Archive and cleanup

**Archive:** `final/archive_manifest.tsv` lists every kept scratch file (19 files, 1.72 GB) with path, bytes, sha256
and role.
- The final tables were checked against `combined_summary.json` before cleanup.
- All 19 files were checked against the manifest again after cleanup, with no mismatch (`final/verification.json`).

**Kept in scratch** (`/global/scratch/users/kh36969/longread_callout/core_context/`):
- `ucd_batch_v0.2/final/`:
  - `all_rows_<kind>.tsv.gz`: every record and reason;
  - `pass_<kind>.tsv.gz`;
  - `strict_lt2500_<kind>.tsv.gz`.
- `ucd_batch_v0.2/remaining/`: catalogue IDs, frozen input table, U/X/D evidence sequences, Check-2 HSPs.
- `ucd_batch_v0.2/rerun/`: input, rows (before/after source), evidence.
- `ucd_batch_v0.1/` (faa63c3): all rows and evidence, unchanged. It is the history, and the source of the reused
  rows' evidence.

**Kept in git:**
- freeze records, source tables (sources.tsv.gz: accession, path, manifest and observed sha256), ID lists;
- the rerun report, summaries, verification and examples;
- the job history and the cleanup lists.

**Removed:** only rebuildable files inside this run's work directories. Each was listed before deletion
(`*/cleanup_list.tsv.gz`, `final/cleanup_extra.json`):
- the rerun's chunk directories and chunk inputs: 4.9 MB;
- the remaining run's chunk directories and chunk inputs (3,800 files, 3.71 GB), rebuildable from
  `input_pairs.tsv.gz`: see `remaining/cleanup.json`;
- the run-level row files of the remaining run, and the `first_batch/` scratch tables (8 files, 0.82 GB). Every row
  is in `final/all_rows_*` (checksummed): see `final/cleanup_extra.json`.

**Total freed:** 4.54 GB.

Shared source FASTA, reads and other tasks' caches were not touched.

## Commands

```bash
D=scale_comamonadaceae_v0.1; S=/global/scratch/users/kh36969/longread_callout/core_context
# rerun of affected pairs (+ controls), from the faa63c3 input table
OUT=$D/ucd_batch_v0.2/rerun WORK=$S/ucd_batch_v0.2/rerun sbatch --export=ALL ucd_batch.sbatch freeze \
    --from-input $D/ucd_batch_v0.1/input_pairs.tsv.gz --ids $D/ucd_batch_v0.2/rerun/rerun_ids.txt --chunk-size 1000
# remaining = boundary catalogue minus processed IDs
OUT=$D/ucd_batch_v0.2/remaining WORK=$S/ucd_batch_v0.2/remaining sbatch --export=ALL --mem=90G ucd_batch.sbatch freeze \
    --catalogue boundary --exclude-input $D/ucd_batch_v0.1/input_pairs.tsv.gz
OUT=... WORK=... sbatch --export=ALL --array=0-759%64 ucd_batch.sbatch chunk    # resubmit to resume
OUT=... WORK=... sbatch --export=ALL --mem=90G ucd_batch.sbatch merge
python ucd_combine.py rerun-report ...     # see rerun/rerun_report.json
sbatch ucd_combine_final.sh                # final combine (ucd_combine.py combine)
python ucd_batch_examples.py --rows $S/ucd_batch_v0.2/final/strict_lt2500_candidate_insertion.tsv.gz ...
```
