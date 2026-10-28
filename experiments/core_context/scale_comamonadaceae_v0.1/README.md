# Assembly-only scale-up v0.1, batch 1 (Comamonadaceae pool)

Code: `../scale_entry.py`, `../scale_compare.py`, `../scale_batch.sbatch`.
Tests: `../test_scale_compare.py`, `../test_scale_entry.py`.
Work files (genes, families, occurrences, compare shards) are on scratch:
`/global/scratch/users/kh36969/longread_callout/core_context/scale_comamonadaceae_v0.1`.

- Evidence is `assembly_only`, `production_ready=false`.
- The entry is the full protein catalog (every GFF3 CDS / FAA protein). Gaps and empty sites are not inputs, and
  reads play no part in selecting assemblies.
- Counts below are first-pass protein-entry candidates (the K=3 neighbour rule), not reviewed U–core–D candidates
  and not independent insertion events.

## Pool and batches (`manifest.tsv`)

- **Pool:** 1,682 verified assemblies from the frozen S2 catalog: FNA + GFF3 + FAA, `verify.json` verified, not
  atypical, every assembly level. 23 assemblies without GFF3/FAA are left out.
- **Hashes:** each file's sha256 is taken from its `verify.json` and re-checked when the genes step reads it.
- **Batches:**

| batch | assembly level | assemblies |
|---:|---|---:|
| 1 | Complete Genome | 412 |
| 2 | Chromosome + Scaffold | 460 |
| 3 | Contig | 810 |

- The 16 subset assemblies are all in batch 1, with identical paths and hashes.

## Design choices

**One protein-family mapping for all batches.** All 4,676,227 distinct proteins of the 1,682 assemblies (from
8,566,868 FAA records) were clustered once at the frozen S4 settings: family 0.9/0.8, homology group 0.5/0.8, with
the member check and repair. This happens before any batch is compared, so a family never depends on batch order.
- Member check: 433 family and 1,741 homology-group members failed against their representative and were repaired
  by the existing rule (2,174 repairs, thresholds unchanged).
- Result: 1,761,494 families and 460,379 homology groups.
- The incremental mode (batch 1 clustered, later batches searched against existing representatives) is implemented
  and tested, but not used. On the smoke run (the 5 comamonas subset assemblies split 3 + 2) it changed 10 of 242
  occurrences with a different-context partner (7 lost, 3 gained) compared with one-pass clustering. That is a
  batch-boundary effect, so it was replaced by the one-pass pool mapping.

**Exact pair accounting with no cap.** `protein_context.compare` kept only the first 20,000 pairs per core in input
order. `scale_compare` replaces it for the scale-up:
- Occurrences collapse to their K=3 context signatures. Every occurrence pair is counted exactly once through
  signature pairs, split into same-assembly and cross-assembly pairs.
- Work runs in resumable blocks. A pair in an unfinished block is `pairs_not_compared`: never similar, different or
  absent.
- Tested:
  - equal to `protein_context.relation` on every pair of random cores;
  - a core above the old cap;
  - one copy per assembly across assemblies;
  - a two-gene core;
  - an interrupted and resumed shard that gives byte-identical output.
- Both structures (single-gene, two-gene) and both scopes (`all`, `target_lt2500`) are compared. Multi-copy
  within an assembly and single-copy across assemblies are both kept. Discovery never requires two copies in one
  assembly.

## Batch 1 results (`summary_b1.json`)

**Completeness:** 0 shards missing, and every core is complete.

| scope : structure | cores (≥ 2 occ.) | pairs total = accounted | not compared | cores above the old 20,000 cap | largest core (occ. / signatures) |
|---|---:|---:|---:|---:|---|
| target_lt2500 : single_gene | 227,694 | 28,600,469 | 0 | 64 | 641 / 173 |
| target_lt2500 : two_gene | 193,481 | 11,996,387 | 0 | 6 | 386 / 87 |
| all : single_gene | 234,342 | 29,502,955 | 0 | 66 | 641 / 216 |
| all : two_gene | 251,662 | 16,039,926 | 0 | 6 | 386 / 87 |

**Candidate counts** (target_lt2500; cores with ≥ 1 `different_context_candidate` pair):

| | single-gene | two-gene |
|---|---:|---:|
| cores with a different-context pair | 5,537 | 1,373 |
| … with a cross-assembly different pair | 5,047 | 1,268 |
| … with a same-assembly different pair | 1,135 | 224 |
| … cross-assembly only | 4,402 | 1,149 |
| … one copy per assembly everywhere | 4,104 | 1,080 |
| different-context occurrence pairs (same-assembly) | 545,804 (15,830) | 46,598 (1,427) |

How to read these counts:
- **Most first-pass candidates are not mobile-element cores.** By product keyword, the 5,537 single-gene cores are:
  other 3,839; hypothetical 980; transposase / IS 527; integrase / recombinase 114; translation / ribosomal 35;
  phage 29; intron RT 13.
- **The largest counts are conserved genes across distant genera.** They include EF-Tu (641 occurrences, 97,966
  different pairs), GroEL (one copy in each of 295 assemblies, 21,184 cross-assembly different pairs) and ETF.
  The K=3 rule cannot separate movement from paralogy or from inheritance with rearrangement. The same-assembly
  split, the S6 site groups and the boundary review are still needed per core.
- **The counts are pairs and cores, not positions or events.**

**Regression on the frozen candidates** (`regression_b1.tsv`; each pair located by coordinates and rechecked at
pool scale):
- **Reproduced:** all 12 main candidates (UCD01–12) and UCDR01, 04, 05 and 06. Each pair is in the same pool core
  and is `different_context_candidate`, and both copies have a different-context partner.
- **UCDR02 (PaaK) changed.** At pool scale its two copies (protein identity 90.3%, coverage 98.9%, at the 0.9 family
  threshold) fall into different families, so it is no longer a family-level core pair. This is a threshold-edge
  effect of clustering 4.7 M proteins, not a missed comparison.
- **UCDR03 (IS30) is unchanged, not lost.** Its K=3 relation is `insufficient_context_labels` here, and it was
  already that in the subset run (`regression_core_check.tsv`). The nearest neighbours of copy B are pseudogene IS5
  fragments without protein labels. It was frozen through the S6 host windows, which skip mobile genes, not
  through the K=3 rule.

**Direct core check of the frozen pairs** (`../subsets_20261027/ucd_candidates_v0.1/regression_core_check.tsv`,
`pair_core_check.py`):
- Main 12: 7 have identical core DNA. Core DNA identity is 94.342–100% and protein identity 97.037–100%, with 100%
  DNA coverage. All 12 are `different_context_candidate` at K=3.
- UCDR03 now has the direct comparison that core_review lacked: protein 99.741%, DNA 99.828%, coverage 100%.

## Resources (batch 1; `time/` on scratch)

| step | tasks | wall time | peak RSS |
|---|---:|---|---|
| pool proteins (all 1,682 FAA) | 1 | 2 min 50 s | 2.6 GB |
| families (pool, 40 CPUs) | 1 | 24 min 22 s | 11.1 GB |
| genes (412 assemblies) | 16 | 30 s max, 393 s total | 0.08 GB |
| occurrences (4.26 M: 2,157,728 single + 2,105,194 two-gene) | 1 | 14 min 2 s | 4.6 GB |
| compare (64 shards, both scopes) | 64 | 23 s max, 1,073 s total | 0.23 GB |
| summarize | 1 | 50 s | 2.7 GB |

Scratch use: families 14 GB, pool proteins 1.1 GB, genes 0.65 GB, occurrences 0.56 GB, compare 0.35 GB.

Provenance notes:
- The pool-protein job started a few minutes before commit `79f8445` and logs `9d825b0` plus uncommitted changes;
  the code it ran is the code committed in `79f8445`.
- The summarize step was rerun after fixing a lookup bug in the regression table (`sig` vs `sig[k]`), which had
  reported every copy as having no different partner. The comparison outputs were not affected.

Not done yet:
- batches 2 and 3;
- S6 site groups and boundary review at scale;
- long-read checks.

Stopped here for review.
