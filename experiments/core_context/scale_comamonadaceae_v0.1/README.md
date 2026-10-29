# Assembly-only scale-up v0.1 (Comamonadaceae pool)

Code: `../scale_entry.py`, `../scale_compare.py`, `../scale_review.py`, `../scale_batch.sbatch`.
Tests: `../test_scale_compare.py`, `../test_scale_entry.py`, `../test_scale_review.py`.
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

## Per-occurrence partner fix (`3634fa5`)

**The bug:** compare v1 kept one cross-assembly and one same-assembly different-context partner per context
signature. A signature merges occurrences from several assemblies, so that single partner did not hold for every
member.

**The fix:** v2 keeps partners per assembly and names them per occurrence (`diff_partners_cross`,
`diff_partners_same_assembly`, `diff_partner_cross`, `diff_partner_same`). It is tested against brute force on every
occurrence and on the review's counterexample.

**Batch-1 impact** (target scope; pair counts and relations unchanged, no v1 partner was wrong):

| missed by v1 | single-gene | two-gene |
|---|---|---|
| named same-assembly partner | 1,197 of 7,965 occurrences (125 cores) | 127 of 1,181 (24 cores) |
| named cross-assembly partner | 100 of 65,298 (79 cores) | 25 of 11,429 (19 cores) |

Batch 1 was recomputed with v2. The batch 1+2+3 joint compare below uses v2.

## Batch 1+2+3 joint compare (`summary_b1-2-3.json`, `regression_b1-2-3.tsv`)

All 1,682 assemblies are compared together under the one pool mapping. This is a joint comparison, not three per-batch
runs: every cross-batch pair is counted.

| target_lt2500 | single-gene | two-gene |
|---|---:|---:|
| cores (≥ 2 occ.) | 729,832 | 671,547 |
| pairs total = accounted | 347,464,545 | 137,843,165 |
| pairs not compared | 0 | 0 |
| cores above the old 20,000 cap | 2,625 | 729 |
| cores with a different-context pair | 17,909 | 4,836 |
| … cross-assembly only | 15,761 | 4,461 |
| … with a same-assembly different pair | 2,148 | 375 |
| different-context pairs (same-assembly) | 5,100,551 (30,169) | 414,711 (4,228) |

- Frozen-candidate regression is the same as batch 1: 12/12 main and UCDR01, 04, 05, 06 reproduced, UCDR02 split into
  two families, UCDR03 `insufficient_context_labels` at K=3.
- These are first-pass K=3 counts. S6 and boundary review have so far run on batch 1 only.

## Batch-1 candidate review: S6, physical overlap, boundaries (`s6_summary_b1.json`)

**Inputs:**
- The S6 review uses the original ordered, oriented host windows: N=5 genes per side, mobile-keyword genes skipped, at
  most 15 scanned. It does not use the compressed K=3 label sets.
- Every target core with a pair that is not `similar_context` is reviewed. That includes cores the K=3 pass left as
  `insufficient_context_labels` or edge cases, such as UCDR03 (IS30).

**Two independent dimensions per pair:**

| dimension | question | field |
|---|---|---|
| pair structure | Are the two copies the same core, with host context different on both sides and no window-overlap artefact? | `ucd_structural_pair` (same assembly: context different on both sides, no host gene shared physically); direct core comparison and boundary review |
| S6 grouping | Which position does each copy correspond to in other assemblies, and can it be counted once? | `s6_grouping` = `supported_positions` (both copies in a site group or single-observation site) or `grouping_unresolved`, with `s6_status_a/b` |

- **Grouping never vetoes structure.** Grouping status depends on which other assemblies are in the core. The pair's
  own structure does not.
- **Version record of the old rule:** strict rule v0 (structure *and* supported positions) gave 12,266 + 825 pairs.
  Its results are kept in `strict_v0_b1/` and as the column `ucd_pair_candidate_strict_v0`.
- **Regression test** (`test_other_occurrences_change_s6_grouping_but_not_the_pair_structure`): adding other
  occurrences turns a pair's grouping from supported to unresolved, while its context class, supports, physical overlap
  and structural call stay the same.
- **No site or event counts from unresolved copies.** Grouping-unresolved copies are never used to count positions or
  independent events.
- **All counts are pairs or occurrences,** not sites or insertion events.

### S6 status (all occurrences of the reviewed cores)

| | single-gene | two-gene |
|---|---:|---:|
| reviewed cores | 40,870 | 21,225 |
| … K=3 had no different pair | 35,333 | 19,852 |
| occurrences | 563,787 | 229,372 |
| site group / single observation / unresolved | 227,048 / 6,634 / 330,105 | 113,236 / 1,635 / 114,501 |
| unresolved: conflict in component | 269,980 | 89,474 |

### Same-assembly U–core–D structural pairs

| | single-gene | two-gene |
|---|---:|---:|
| same-assembly pairs | 46,408 | 6,851 |
| context different on both sides | 20,616 | 2,068 |
| … minus pairs sharing host genes physically | −128 | −19 |
| **structural pairs** | **20,488** | **2,049** |
| … S6 supported positions (= strict v0) | 12,266 | 825 |
| … S6 grouping unresolved | 8,222 | 1,224 |
| cores with a structural pair | 1,259 | 290 |
| … only grouping-unresolved pairs | 369 | 81 |
| … K=3 had no different pair | 257 | 94 |
| cores by product: transposase/IS · other · hypothetical · integrase/recombinase · intron RT · phage · translation | 642 · 406 · 163 · 29 · 12 · 1 · 6 | 126 · 114 · 40 · 8 · 0 · 2 · 0 |

The physical-overlap pairs are not dropped: they stay in `s6_cores_b1.tsv.gz` with their context class. Single-gene
pairs sharing host genes physically: 306 context repeated both sides, 128 context different both sides, 24 repeated
one side, 682 context unknown.

### Direct core comparison and boundaries, same-assembly structural pairs

Pair columns are in `boundary_pairs_same_assembly_b1.tsv.gz`; occurrence columns in
`boundary_occurrences_same_assembly_b1.tsv.gz`.

| | single: supported | single: unresolved | two-gene: supported | two-gene: unresolved |
|---|---:|---:|---:|---:|
| pairs | 12,266 | 8,222 | 825 | 1,224 |
| identical core DNA | 5,447 | 5,392 | 444 | 823 |
| core DNA identity min | 75.373 | 78.553 | 84.345 | 86.278 |
| core DNA coverage < 0.9 | 251 | 55 | 16 | 16 |
| no core DNA hit | 0 | 0 | 0 | 0 |
| both core ends covered | 11,649 | 7,613 | 689 | 1,076 |

- **Low-identity pairs are not hidden.** Identity and coverage are reported and ranked, never used to admit or drop a
  pair. Pairs below 90% core DNA identity or with coverage < 0.9 are listed.
- **Strict pairs are unchanged:** all 13,091 strict-v0 pairs give identical core and boundary values in this run (one
  all-vs-all per core in place of one per assembly; blast runs with `-subject`, so each pair is scored on its own).
- **Per-occurrence boundary** over all structural same-assembly partners: the boundary summary was recomputed for every
  occurrence, including those whose partner set grew.

| occurrences | single-gene | two-gene |
|---|---:|---:|
| with a structural partner | 9,328 | 1,531 |
| both sides resolved | 6,600 | 1,002 |
| estimated span < 2,500 bp | 5,730 | 501 |
| estimated span ≥ 2,500 bp (listed separately) | 870 | 501 |
| boundary unknown (kept) | 2,728 | 529 |
| … of which with an S6-supported partner: < 2,500 / ≥ 2,500 / unknown | 3,046 / 422 / 1,541 | 317 / 258 / 233 |

### Cross-assembly single-copy pairs

- **Which pairs:** two assemblies, each holding one copy of the core, with host context different on both sides
  (`different_position_candidate`). Neither multi-copy assemblies nor empty sites are required.
- **Same review:** these pairs get the same direct core comparison and boundary review, and the same S6 grouping field.
- **Physical overlap** does not apply across assemblies.
- **Files:** pair rows are on scratch (`review/s6/b1/*/cross_assembly_single_copy_pairs.tsv.gz`,
  `review/boundary/b1/*/pairs.tsv.gz`); per-occurrence boundaries are in
  `boundary_occurrences_cross_single_copy_b1.tsv.gz`.

| | single-gene | two-gene |
|---|---:|---:|
| cross-assembly different-position pairs (all) | 441,805 | 33,323 |
| … both assemblies single-copy | 187,956 | 18,118 |
| … S6 supported / grouping unresolved | 43,478 / 144,478 | 6,809 / 11,309 |
| cores | 3,451 | 753 |
| … only grouping-unresolved pairs | 1,655 | 382 |
| … K=3 had no different pair | 305 | 86 |
| cores by product: other · hypothetical · transposase/IS · integrase/recombinase · translation · phage · intron RT | 2,470 · 660 · 184 · 91 · 28 · 16 · 2 | 591 · 110 · 28 · 12 · 8 · 4 · 0 |
| identical core DNA | 142 | 59 |
| core DNA coverage < 0.9 | 4,078 | 2,197 |
| occurrences: < 2,500 / ≥ 2,500 / unknown | 11,069 / 615 / 26,516 | 2,077 / 413 / 3,540 |

How to read the cross-assembly results:
- **Not a mobile-element set.** Unlike the same-assembly set, few of these pairs have identical core DNA, and most
  cores are not mobile-keyword cores.
- **Different context is ambiguous here.** Across assemblies it also fits inherited genes in rearranged or divergent
  genomes, so the same caution as for the K=3 counts applies.
- **Many boundaries are unknown.** One core end is often `not_reached`, meaning the chain of aligned blocks stops
  before that end in one of the two copies.

### Frozen regression set (`s6_regression_b1.tsv`)

**Structural pairs in both tiers:** all 12 main candidates and UCDR01, 03, 04, 05, 06. Their direct core values equal
`../subsets_20261027/ucd_candidates_v0.1/regression_core_check.tsv`.

| candidate | S6 grouping | copy status (a / b) | core DNA % | protein % |
|---|---|---|---:|---:|
| UCD01 | grouping_unresolved | unresolved_pairs / unresolved_pairs | 100.0 | 100.0 |
| UCD03 | grouping_unresolved | unresolved_pairs / unresolved_pairs | 100.0 | 100.0 |
| UCD07 | grouping_unresolved | unresolved_pairs / unresolved_pairs | 94.342 | 97.037 |
| UCDR01 | grouping_unresolved | conflict_in_component / site_group | 99.412 | 99.242 |
| UCDR04 | grouping_unresolved | conflict_in_component / site_group | 99.328 | 100.0 |
| UCDR06 | grouping_unresolved | site_group / conflict_in_component | 100.0 (cov. 0.9966) | 100.0 |
| other 9 main, UCDR03, UCDR05 | supported_positions | site group or single observation | 99.895–100 | 99.741–100 |

**What changed for the six:**
- In the 16-assembly subset these copies had supported positions.
- At 412 assemblies their copies join cores with many more members. In those cores some pairwise verdicts stay
  unresolved or the site-group component is not a clique.
- Their sequences, families, host contexts and physical-overlap status are unchanged. They remain local U–core–D
  candidates whose cross-sample position grouping is unresolved.

**UCDR03 (IS30):** `insufficient_context_labels` at K=3. In the host-window review it is a structural pair with
supported positions (protein 99.741%, DNA 99.828%).

**UCDR02 (PaaK):** not a core pair at pool scale (split into two families, see above).

## Resources (batch-1 review)

| step | tasks | wall time max |
|---|---:|---|
| windows | 1 | 28 min |
| s6 | 64 | 2 min 40 s |
| boundary (structural + cross pairs) | 64 | 42 min 40 s |
| summarize | 1 | 1 min 34 s |

The first s6 submission failed on an argument-parser bug (`--repo-out` required for s6), fixed in `4175584`.

Not done yet:
- S6 and boundary review on the batch 1+2+3 joint set (windows for batches 2 and 3 are built);
- long-read checks, batch by batch after candidate review; no FASTQ downloaded.

Stopped here for review.
