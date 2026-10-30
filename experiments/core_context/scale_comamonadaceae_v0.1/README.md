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

### Definitions used in the evidence tables

Three layers are reported separately, never merged into one "candidate" label:
1. **Context different:** the pair entered structural review.
2. **Core match:** how well the two cores align.
3. **Length target:** the whole element is resolved on both sides with an estimated span under 2,500 bp. A short coding
   core does not stand in for this.

**Core match class** (`core_match_class`; reporting and ranking only, not an admission rule):

| class | definition |
|---|---|
| `identical_core_dna` | the two coding spans (core-oriented DNA) are identical strings |
| `homolog_full_coverage` | not identical. Best blastn HSP of core vs core (`-task blastn`, e-value 1e-5) covers ≥ 90% of each copy's core DNA; `core_dna_coverage_min` = min(HSP length / core length) over the two copies. No identity threshold beyond family membership (protein 0.9 identity / 0.8 coverage); `core_dna_identity` (pident of that HSP) is reported as a separate column |
| `partial_core_match` | best HSP covers < 90% of at least one copy: a local match, kept for review |
| `no_core_dna_hit` | no core-vs-core blastn HSP |

**Both core ends covered vs boundary resolved:**
- **Both core ends covered** is a pair property: the chained ±3 kb flank alignment of the two copies covers the upstream
  and the downstream core end in both copies (`covered` or `covered_at_chain_end`). Only then is the shared extension
  past each end measurable for that pair; otherwise it is unknown, never 0.
- **Boundary resolved** is an occurrence property within one pair scope:
  - every partner of that scope gives a measured extension on that side;
  - the extensions agree within 20 bp;
  - the shared stretch does not reach the window or contig end.

  When both sides are resolved, estimated span = coding span + median upstream + median downstream extension, then
  classed `lt2500` or `ge2500`. Anything else is `unknown`.

**Counting scopes:**
- Pairs, distinct occurrences, cores and S6 status are counted separately.
- An occurrence can be in both the same-assembly and the cross-assembly scope. Its boundary is computed per scope, so
  occurrence counts of the two scopes are not added.

### Cross-assembly pairs (every copy class)

- **Which pairs:** every pair of copies from two assemblies with host context different on both sides
  (`different_position_candidate`), whatever the copy numbers.
- **Copy class** (single×single, single×multi, multi×multi) is a classification only.
- **Why multi-copy pairs are needed:** two copies with similar contexts in one assembly give no same-assembly pair,
  yet each can differ from a copy elsewhere (`test_two_similar_copies_in_one_assembly_differ_from_a_single_copy_elsewhere`).
  Commit `d662ad7` reviewed only single×single; `b7d40ef` adds the rest.
- **Reuse:** single×single rows were reused from that run (`--reuse`). Reuse is exact: a recomputed shard matched row
  for row, because blast scores each pair on its own.
- **Physical overlap** does not apply across assemblies.
- **Files:**
  - pair rows are on scratch (`review/boundary/b1/*/pairs.tsv.gz`);
  - per-occurrence boundaries: `boundary_occurrences_cross_assembly_b1.tsv.gz`.

| single-gene | single×single | single×multi | multi×multi | all |
|---|---:|---:|---:|---:|
| pairs | 187,956 | 54,358 | 199,491 | 441,805 |
| … S6 supported / grouping unresolved | 43,478 / 144,478 | 13,523 / 40,835 | 54,692 / 144,799 | 111,693 / 330,112 |
| cores (a core can be in several classes) | 3,451 | 797 | 526 | 4,047 |
| identical core DNA | 142 | 813 | 23,007 | 23,962 |
| full-coverage homolog | 183,734 | 52,599 | 175,596 | 411,929 |
| partial core match | 4,078 | 946 | 888 | 5,912 |
| no core DNA hit | 2 | 0 | 0 | 2 |
| both core ends covered | 100,387 | 33,693 | 137,513 | 271,593 |
| core DNA identity min | 73.063 | 74.71 | 75.653 | 73.063 |

| two-gene | single×single | single×multi | multi×multi | all |
|---|---:|---:|---:|---:|
| pairs | 18,118 | 4,315 | 10,890 | 33,323 |
| … S6 supported / grouping unresolved | 6,809 / 11,309 | 1,137 / 3,178 | 1,686 / 9,204 | 9,632 / 23,691 |
| cores | 753 | 165 | 114 | 912 |
| identical core DNA | 59 | 51 | 437 | 547 |
| full-coverage homolog | 15,862 | 3,784 | 8,689 | 28,335 |
| partial core match | 2,197 | 480 | 1,764 | 4,441 |
| both core ends covered | 9,680 | 1,512 | 4,710 | 15,902 |

How to read the cross-assembly results:
- **Multi-copy pairs carry most identical cores.** Identical core DNA rises from 142 (single×single) to 23,962
  single-gene pairs. A multi-copy core gives many pairs, so this is not the number of new sites or families.
- **Different context is ambiguous across assemblies.** It also fits inherited genes in rearranged or divergent
  genomes.

**Per-occurrence boundaries** (cross-assembly scope, over all cross-assembly partners of the occurrence):

| | single-gene | two-gene |
|---|---:|---:|
| occurrences | 50,562 | 7,774 |
| both sides resolved | 15,719 | 3,110 |
| estimated span < 2,500 / ≥ 2,500 / unknown | 14,544 / 1,175 / 34,843 | 2,350 / 760 / 4,664 |

**Distinct counts across the two scopes** (batch 1):

| | single-gene | two-gene |
|---|---:|---:|
| occurrences: same-assembly / cross-assembly / in both / distinct | 9,328 / 50,562 / 7,275 / 52,615 | 1,531 / 7,774 / 1,131 / 8,174 |
| cores: same-assembly / cross-assembly / distinct | 1,259 / 4,047 / 4,573 | 290 / 912 / 1,045 |

Same-assembly core match classes (all structural pairs):
- single-gene: 10,839 identical, 9,343 full-coverage homolog, 306 partial;
- two-gene: 1,267 identical, 750 full-coverage homolog, 32 partial.

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
| boundary (structural + single×single cross pairs) | 64 | 42 min 40 s |
| s6 rerun (all cross-assembly pairs) | 64 | 25 s |
| boundary rerun (adds single×multi and multi×multi, reuses the rest) | 64 | 24 min 37 s |
| summarize | 1 | 1 min 34 s |

The first s6 submission failed on an argument-parser bug (`--repo-out` required for s6), fixed in `4175584`.

## Batch 1+2+3 joint review (`s6_summary_b1-2-3.json`, `candidates_summary_b1-2-3.json`)

The rules, thresholds and definitions are the same as for batch 1 (see "Definitions used in the evidence tables"). The
results are `assembly_only`, `production_ready=false`. Pairs, occurrences and cores are not sites or insertion events.

**Completeness:**
- S6, boundary and candidates ran over all 64 shards; `shards_missing` is empty.
- `pairs_pending_review` is 0: every same-assembly structural pair and every cross-assembly different-context pair was
  reviewed directly (core and boundary).
- Batch-1 rows were reused (`REUSE=b1`): 20,488 + 441,805 single-gene and 2,049 + 33,323 two-gene pairs.

### S6 status and pair tiers

| | single-gene | two-gene |
|---|---:|---:|
| reviewed cores | 203,268 | 139,616 |
| occurrences: site group / single observation / unresolved | 1,653,640 / 9,756 / 2,927,983 | 1,168,654 / 2,788 / 1,372,537 |
| same-assembly structural pairs | 37,898 | 4,865 |
| … S6 supported / grouping unresolved | 18,189 / 19,709 | 2,187 / 2,678 |
| … cores (with an S6-supported pair) | 2,188 (1,160) | 455 (260) |
| cross-assembly different-context pairs | 3,579,311 | 272,843 |
| … S6 supported / grouping unresolved | 490,152 / 3,089,159 | 56,605 / 216,238 |
| … cores (with an S6-supported pair) | 13,149 (6,226) | 3,128 (1,393) |

### Cross-assembly pairs by copy class

| single-gene | single×single | single×multi | multi×multi | all |
|---|---:|---:|---:|---:|
| pairs | 2,358,681 | 440,031 | 780,599 | 3,579,311 |
| … S6 supported / grouping unresolved | 278,088 / 2,080,593 | 55,933 / 384,098 | 156,131 / 624,468 | 490,152 / 3,089,159 |
| cores (a core can be in several classes) | 12,012 | 2,269 | 1,165 | 13,149 |
| identical core DNA | 677 | 4,201 | 55,469 | 60,347 |
| full-coverage homolog | 2,307,285 | 426,237 | 708,591 | 3,442,113 |
| partial core match | 50,609 | 9,593 | 16,539 | 76,741 |
| no core DNA hit | 110 | 0 | 0 | 110 |
| both core ends covered | 1,193,837 | 269,854 | 497,679 | 1,961,370 |
| core DNA identity min | 70.068 | 73.734 | 75.653 | 70.068 |

| two-gene | single×single | single×multi | multi×multi | all |
|---|---:|---:|---:|---:|
| pairs | 183,483 | 39,625 | 49,735 | 272,843 |
| … S6 supported / grouping unresolved | 41,585 / 141,898 | 5,613 / 34,012 | 9,407 / 40,328 | 56,605 / 216,238 |
| cores | 2,889 | 413 | 218 | 3,128 |
| identical core DNA | 312 | 171 | 7,524 | 8,007 |
| full-coverage homolog | 157,081 | 34,884 | 34,723 | 226,688 |
| partial core match | 26,090 | 4,570 | 7,488 | 38,148 |
| both core ends covered | 104,877 | 11,457 | 18,414 | 134,748 |
| core DNA identity min | 74.597 | 75.163 | 78.649 | 74.597 |

**Same-assembly core match classes** (structural pairs):

| | identical | full-coverage homolog | partial | both core ends covered |
|---|---:|---:|---:|---:|
| single-gene, S6 supported | 6,735 | 10,666 | 788 | 16,647 |
| single-gene, grouping unresolved | 10,454 | 8,857 | 398 | 17,831 |
| two-gene, S6 supported | 1,539 | 620 | 28 | 2,009 |
| two-gene, grouping unresolved | 1,748 | 882 | 48 | 2,238 |

Full-coverage homologs go down to 70–76% core DNA identity. Coverage, not identity, is the criterion, so the identity
column must be read alongside the class.

### Per-occurrence boundaries and distinct counts

"Both core ends covered" counts pairs. "Both sides resolved" counts occurrences, per scope, over all partners in that
scope.

| | single-gene same | single-gene cross | two-gene same | two-gene cross |
|---|---:|---:|---:|---:|
| occurrences | 18,486 | 259,174 | 3,041 | 37,956 |
| both sides resolved | 12,237 | 65,904 | 1,813 | 14,298 |
| span < 2,500 / ≥ 2,500 / unknown | 10,745 / 1,492 / 6,249 | 61,759 / 4,145 / 193,270 | 938 / 875 / 1,228 | 11,172 / 3,126 / 23,658 |
| with an S6-supported partner | 7,637 | 63,149 | 1,352 | 11,737 |

| distinct across scopes | single-gene | two-gene |
|---|---:|---:|
| occurrences: same / cross / in both / distinct | 18,486 / 259,174 / 16,286 / 261,374 | 3,041 / 37,956 / 2,563 / 38,434 |
| cores: same / cross / in both / distinct | 2,188 / 13,149 / 1,581 / 13,756 | 455 / 3,128 / 295 / 3,288 |

### Candidate cores: matched core, different context, both spans < 2,500 bp

**Script:** `scale_candidates.py` (stage `candidates`).

**Pair evidence state:** each structural or cross-assembly pair gets one state, in the scope of that pair:

| state | meaning |
|---|---|
| `match_both_lt2500` | core match is `identical_core_dna` or `homolog_full_coverage`, and both copies have span class `lt2500` |
| `match_span_ge2500` | sufficient core match, at least one copy `ge2500` |
| `match_boundary_unknown_pending` | sufficient core match, at least one copy `unknown`, none `ge2500` |
| `partial_core_match_pending` | best core HSP covers < 90% of a copy |
| `no_core_dna_hit` | no core-vs-core HSP |

Pending pairs are kept and counted per core. They are not dropped and not promoted.

| pairs | single-gene same | single-gene cross | two-gene same | two-gene cross |
|---|---:|---:|---:|---:|
| `match_both_lt2500` | 16,177 | 70,400 | 1,069 | 7,682 |
| `match_span_ge2500` | 4,747 | 12,280 | 2,446 | 13,317 |
| `match_boundary_unknown_pending` | 15,788 | 3,419,780 | 1,274 | 213,696 |
| `partial_core_match_pending` | 1,186 | 76,741 | 76 | 38,148 |
| `no_core_dna_hit` | 0 | 110 | 0 | 0 |

**Core list** (`core_candidates_b1-2-3.tsv.gz`): one row per core with at least one `match_both_lt2500` pair, listing:
- its pair counts by scope and copy class;
- S6-supported pairs, identical pairs and identity min;
- occurrence S6 status;
- pending counts;
- best pair;
- local-read assemblies.

The pairs themselves are in `match_both_lt2500_pairs_b1-2-3.tsv.gz`.

| | single-gene | two-gene |
|---|---:|---:|
| cores | 5,416 | 1,245 |
| … with a same-assembly pair / cross-assembly only | 1,346 / 4,070 | 221 / 1,024 |
| … with an S6-supported qualifying pair | 2,588 | 590 |
| … K=3 had no different pair for this core | 591 | 172 |
| … transposase/IS, integrase/recombinase/resolvase | 830, 88 | 91, 11 |
| … hypothetical, other | 1,472, 2,969 | 302, 831 |
| distinct occurrences in a qualifying pair | 32,629 | 6,118 |
| distinct occurrences only in pending pairs | 227,141 | 30,743 |

**Housekeeping paralogs rank high.** The sort is: identical core first, then local reads, then pair count. This puts
duplicated housekeeping genes near the top:
- elongation factor Tu (`tufA`/`tufB`, 570 same-assembly pairs);
- Hcp (T6SS).

They have the same core and different flanks by construction. Product is for reporting and ranking only, so they stay
in the list; reviewers should read `product_hint`.

### Validation batch 1 (`validation_batch1_b1-2-3.tsv.gz`)

**What it contains:** occurrences in `match_both_lt2500` pairs whose assembly has verified local long reads. Verified
means `data/reads/<run>/verify.json` has status verified and an md5 match. Runs are linked through `runs.tsv`; 151
assemblies have them.

**Columns:** coordinates, estimated span per scope, partner counts, whether a partner assembly also has reads, and runs
with their platforms.

Nothing was downloaded and no reads were mapped.

| | single-gene | two-gene |
|---|---:|---:|
| occurrences | 5,756 | 1,000 |
| cores | 1,759 | 437 |
| assemblies | 149 | 141 |

**Suggested order:** ad hoc tiers, not written to the file.
- **A:** the occurrence has a qualifying same-assembly partner, so one read set covers both copies.
- **B:** cross-assembly only, and a partner assembly also has reads.
- **C:** cross-assembly only, and only this side has reads.

"Identical" below means the core has at least one identical-core pair.

| tier | single-gene occ. / cores / assemblies | two-gene occ. / cores / assemblies |
|---|---|---|
| A, identical | 1,661 / 226 / 141 | 63 / 20 / 15 |
| A, homolog only | 707 / 166 / 113 | 88 / 38 / 30 |
| B, identical | 162 / 72 / 75 | 24 / 9 / 17 |
| B, homolog only | 1,122 / 274 / 131 | 260 / 86 / 102 |
| C, identical | 57 / 44 / 34 | 10 / 6 / 10 |
| C, homolog only | 2,047 / 1,068 / 149 | 555 / 294 / 127 |

**Platforms** (occurrence × run, single-gene):
- ONT: 4,502;
- PacBio CLR subreads: 391;
- PacBio HiFi/CCS: 257;
- PacBio with read type unresolved: 1,243.

### Frozen regression set (`s6_regression_b1-2-3.tsv`)

**Structural pairs:** all 12 main candidates and UCDR01, 03, 04, 05, 06 are still structural pairs. UCDR02 is still not a
core pair.

**S6 grouping:**
- UCD02, 04, 05, 06, 08, 09, UCDR03 and UCDR05: `supported_positions`.
- UCD01, 03, 07, UCDR01, 04 and 06: `grouping_unresolved`, as in batch 1.
- UCD10, 11 and 12 moved from `supported_positions` (batch 1) to `grouping_unresolved`. One copy of each is now
  `conflict_in_component`.

This is the case covered by `test_other_occurrences_change_s6_grouping_but_not_structure`: more assemblies change the
grouping, not the pair.

### Resources (joint review)

| step | tasks | wall time max |
|---|---:|---|
| s6 | 64 | 5 min 56 s |
| boundary (reuses batch-1 rows) | 64 | 3 h 42 min |
| summarize | 1 | 5 min 36 s |
| candidates | 1 | 2 min 3 s |

Not done yet:
- `match_boundary_unknown_pending` and `partial_core_match_pending` pairs stay pending;
- no FASTQ downloaded.

## Formal validation list (`validation_pairs_b1-2-3.tsv.gz`, `ucd_readcheck.py pairs`)

**What is in it:** one row per `match_both_lt2500` pair with at least one side in an assembly that has verified local
reads (usable run, `verify.json` verified, MD5 equal to the archive). Each row fixes:
- both occurrences and their coordinates;
- the element interval per side: coding span plus that pair's own shared extensions (`a_up`/`a_down`), in genome
  coordinates;
- core match class, core DNA identity and coverage;
- both sides' boundary status in the pair's scope;
- the runs.

**Tiers** are per pair, not per core family:
- **A:** same assembly, so one read set covers both positions.
- **B:** cross assembly, with reads on both sides.
- **C:** cross assembly, with reads on one side only, so at most single-side support.

**Ranking within a tier:** identical core first, then core DNA identity, then coverage.

| pairs | A | B | C |
|---|---:|---:|---:|
| single-gene | 3,151 | 2,732 | 19,240 |
| two-gene | 137 | 247 | 1,837 |
| assemblies with reads | 142 | 138 | 149 |

**Cross-check:** the occurrences with reads in this list are exactly the 6,756 of `validation_batch1_b1-2-3.tsv.gz`.

**Controls** stay in the list as labels, not exclusions:
- **`control_class=conserved_paralog_control`:** product hint `translation_or_ribosomal`, for example elongation factor
  Tu and GreA/GreB. Pairs: A 123, B 2, C 22.
- **`validation_controls_b1-2-3.tsv.gz`:** same-assembly pairs with `context_repeated_both_sides` in assemblies with
  reads (1,499 single-gene, 588 two-gene). These are not U-core-D pairs and are tested on the coding span.

## Long-read check (`readcheck_round1/`, `readcheck_round2/`)

**Question:** do the assembly's own long reads support the local U-core-D structure it shows at each tested
interval?
- The reads are mostly the ones the assembly was built from, so this is a self-consistency check under the current
  alignment, locus-competition and counting rules.
- It does not exclude every misassembly and does not fix exact biological insertion boundaries.
- It says nothing about mobility or independent insertion events.

**Method (`ucd_readcheck.py`, s3_readcheck v0.3 counting unchanged):**
- `genome`: all reads of every verified run are mapped once to the whole assembly (secondaries kept; SHA256 rechecked
  against `verify.json`). The PAF is kept on scratch (`core_context/readmap/<sample>/`) and reused by every round.
- `check`: reads with a genome alignment within 2 kb of a tested interval are taken from that full-genome PAF. Two
  alleles are built from the same assembly:
  - P = U + element + D, as assembled;
  - A = U + D, element removed.

  No empty site in another assembly is needed.
- **Counting:**
  - a read must span the element plus 300 bp anchors and pass the CIGAR limit;
  - the allele score margin is 20;
  - both anchors must be placed at this position genome-wide; otherwise the read is `off_locus` or
    `locus_ambiguous`;
  - molecule grouping and length-matched counts apply;
  - verdict: insufficient if fewer than 3, consistent if the other allele is at most 10%, otherwise mixed or
    contradicted.
- **Interval class:** supported, conflict, locus_ambiguous, insufficient_coverage, not_tested (interval at a contig
  end) or not_checked (no result).
- **Completeness:**
  - `status` checks every assembly of a round against its expected intervals. `status.tsv` keeps every SLURM attempt
    with its node and error.
  - `summarize` refuses a round with missing results; with `--allow-incomplete`, missing intervals stay `not_checked`
    and outside every count.
  - A failed task is never read as missing support.
- **No downloads:** this check has no fetch stage. It reads only local FASTQ whose `verify.json` and archive MD5
  match.

**Three units, kept apart:**
- **Pair:** a row of the validation list.
- **Occurrence:** an occurrence id.
- **Tested interval:** sample, contig and element interval.
- **How they relate:** one occurrence gets one interval per distinct element boundary its partners give, and one
  interval can carry several occurrence ids. `interval_map.tsv.gz` keeps every (interval, occurrence, pair, side)
  link. `occurrences.tsv.gz` gives each occurrence's intervals and class; `pairs.tsv.gz` gives both sides' interval
  ids and classes.

**Count mode, per interval:**
- `molecules`: one record per molecule, or ZMW from the read name;
- `assumed_molecules`: SRA spot taken as a ZMW;
- `read_records`: no molecule grouping available;
- `mixed`: more than one of these.

Support is reported per mode, never as one depth, and `P_matched_<mode>` / `A_matched_<mode>` are kept per interval.

### Round 1: 9 assemblies, 12 runs

The assemblies were chosen for specific A pairs in both structures and for 3 B assembly pairs. After the unit fix, the
check was rerun from the kept mappings: every count and verdict of the 2,556 intervals is identical to the first run.

| unit | count |
|---|---:|
| pairs | 5,834 |
| occurrences | 1,109 (293 with more than one interval, max 26) |
| tested intervals | 2,556 (32 carry more than one occurrence id) |

All 2,556 intervals are `consistent` and all 1,109 occurrences have every interval supported.

| count mode | intervals | assemblies | P_matched min / median | intervals with any A_matched |
|---|---:|---:|---:|---:|
| molecules | 1,048 | 4 | 9 / 30 | 0 |
| assumed_molecules | 340 | 2 | 16 / 120 | 11 |
| read_records | 1,168 | 3 | 41 / 96 | 0 |

### Round 2: cumulative, all 149 assemblies with local reads and listed intervals, 164 runs

**Selection:** round 1 plus the remaining 140 assemblies. A greedy order by new distinct cores, not pair counts, kept
gaining cores up to about 114 assemblies, and each assembly costs minutes. So all were run, which brings every B pair
to both sides. Round 1's per-assembly results are linked in, after checking that the interval sets are identical.

**Completeness:** 149 / 149 assemblies, 16,641 / 16,641 intervals with results. Four attempts failed with
shared-storage I/O errors and all four were rerun to full results. Whether the cause is the nodes or a short storage
event is open: n0382 failed three times, but the three round-2 failures fell within 20 s and n0334 completed another
task.

| assembly | failed attempt | rerun | state |
|---|---|---|---|
| SAMN08775340 (round 1) | 27441044_3, n0382, `OSError` I/O | 27441820_3 | 435 / 435 |
| SAMN51281044 | 27441725_133, n0334, `OSError` I/O | 27442220_133 | 91 / 91 |
| SAMN51281056 | 27441725_138, n0382, `BrokenPipeError` | 27442220_138 | 61 / 61 |
| SAMN56368954 | 27441725_147, n0382, gzip write | 27442220_147 | 119 / 119 |

| unit | count |
|---|---:|
| pairs | 29,431 |
| occurrences | 9,809 (2,256 with more than one interval, max 26) |
| tested intervals | 16,641 (252 carry more than one occurrence id) |

| interval class | intervals |
|---|---:|
| supported | 16,413 |
| insufficient_coverage | 159 |
| locus_ambiguous | 54 |
| conflict | 2 |
| not_tested (contig end) | 13 |

- **Occurrences:**
  - 9,629 have every interval supported;
  - 178 have no interval supported;
  - 2 have a conflict.
- **Conflict:** both conflicts are the two sides of one repeated-context control, `VC_aa36fecdda17`, in
  SAMEA120399594. It is a 2.4 kb two-gene repeat:
  - one side has 92 P / 34 A (mixed), the other 0 P / 21 A (contradicted);
  - each side has 268–306 `locus_ambiguous` molecules.

  This is consistent with a repeat that the assembly or the reads do not resolve there.
- **Insufficient intervals:**
  - 145 of the 159 `insufficient_coverage` intervals come from two ONT runs with only 1,066 and 1,616 reads mapped to
    the whole assembly (SAMEA115676864, SAMEA115676832): too little coverage, not a failure;
  - `locus_ambiguous` intervals are mainly in SAMEA120399595 / SAMEA120399592, where molecules near the interval
    also place elsewhere.

| count mode, supported intervals | intervals | assemblies | P_matched min / median |
|---|---:|---:|---:|
| molecules | 12,028 | 119 | 3 / 61 |
| assumed_molecules | 1,097 | 11 | 10 / 112 |
| read_records | 3,288 | 19 | 4 / 89 |

Intervals with any A_matched read: 149 molecules, 13 assumed_molecules, 7 read_records.

| pairs | both supported | one supported | none supported | single side supported | single side not supported |
|---|---:|---:|---:|---:|---:|
| A single-gene | 2,996 | 21 | 11 | | |
| A two-gene | 137 | | | | |
| B single-gene | 2,702 | 30 | | | |
| B two-gene | 236 | 9 | | | |
| C single-gene | | | | 19,057 | 165 |
| C two-gene | | | | 1,821 | 12 |
| conserved paralog control | A 122 / B 2 | | A 1 | C 22 | |
| repeated-context control | 1,950 | 130 | 6 | | |

The repeated-context controls also include 1 conflict pair.

**Distinct cores (candidates, controls excluded):** 2,184 cores have a side with reads, and 780 of them (642
single-gene, 138 two-gene) have at least one A or B pair with both intervals supported. The pair totals above are not
counts of independent families.

**Reading:**
- The local assembled structure is supported at 16,413 of 16,641 tested intervals under the current rules.
- Candidates and controls behave alike, apart from the one repeated-context control in conflict. The check therefore
  does not separate a mobile element from a conserved paralog or a segmental duplicate.
- Interpreting the candidates (conserved paralog position, large duplication, different host context, unexplained)
  is a separate layer and has not been done.

**Figures (`readcheck_round*/figures/`, SVG):** one panel per interval, showing U / element / D and the genome
alignments of reads with a call there, colored by call: blue P, red A, orange locus_ambiguous, light grey off_locus,
purple cigar_fail. Pairs shown:
- UCD06 `VP_80c3de3bb257`, UCDR03 `VP_d69542349125`, UCDR04 `VP_cbb7fe69d4b5` (EF-Tu);
- A single-gene `VP_000086f80190`, A two-gene `VP_8c6421b64450`;
- B single-gene `VP_023b5b6ef0fb`, B two-gene `VP_4204ef205cc5`;
- C `VP_0025bbc754af`;
- repeated-context controls `VC_02e5dfdac4fd`, `VC_0081429ee7c1`.

**Files per round:**
- `positions.tsv.gz`: per interval;
- `occurrences.tsv.gz`, `interval_map.tsv.gz`, `pairs.tsv.gz`;
- `runs.tsv`, `status.tsv`, `selection.json`, `samples.txt`, `summary.json`.

Not done yet:
- candidate interpretation layer;
- old pipelines that still fetch inside lr6 compute jobs (`batch_sample.sbatch`, `s3_readcheck.sbatch`,
  `read_callout/01_fetch.sbatch`, `pacbio_hifi_pilot/pilot.sbatch`) are not changed;
- no new downloads.

Stopped here for review.
