# DNA-first U-core-D entry: validation on 50 assemblies (2026-11-09)

DNA entry `dna_discovery.py` (v0.1.2 and v0.2.x); two checks unchanged: `ucd_compare.py` (sha 4909ad16) and
`ucd_checks.py` (sha 35111726), DECISION_VERSION decision_v0.2. Strict rule: pass both checks and the extracted X of
both sides 0 < length < 2500 bp. No protein, gene, empty-site, gap, S6 or read condition.

Assemblies: `validation_manifest.tsv`. 50 Comamonadaceae (batch 1): 40 chosen as dense in protein-entry strict pairs,
10 random complete. All sha256-verified FNA.

## Versions

| version | change | why |
|---|---|---|
| v0.1 | first run | identity filtered on approximate PAF counts (superseded) |
| v0.1.1 | minimap2 `-c` | base-level identity |
| v0.1.2 | occurrence search `-P` (no cap); containment clustering of library | `-N 1000` capped secondary chains per query contig and silently dropped copies of high-copy segments; linclust left the library 3x redundant |
| v0.2.0 | `enumerate` / `align` replace `pairs` | one alignment per canonical key (grid windows, relative orientation, A-B = B-A), cached; one subject per blastn call (results independent of batch); all pairs of high-copy segments enumerated in 20,000-pair chunks, deferred chunks resumable |
| v0.2.1 | keys sharded by subject window | 4.0 M blastn calls -> 0.22 M; per-key results unchanged |

## Results (50 assemblies)

| | v0.1.2 | v0.2 |
|---|---|---|
| source pairs computed | 7,211,086 | 10,493,703 |
| high-copy pairs deferred | 3,282,617 | 0 |
| distinct resolved pairs checked | 604,941 | 854,446 (565,908 reused by fingerprint, 288,538 checked) |
| source check pass | 604,941 / 604,941 | 854,446 / 854,446 |
| both checks pass | 114,850 | 348,575 |
| strict < 2500 pairs (same / cross assembly) | 83,912 (16,727 / 67,185) | 313,666 (23,585 / 290,081) |
| strict distinct X intervals / assemblies | 30,026 / 50 | 42,094 / 50 |
| pair alignment cost (task-hours) | 31.4 (8 CPU tasks) | v0.2.0 62.1, v0.2.1 2.46 (16 CPU tasks) |

Regression against the protein entry's strict pairs between the same 50 assemblies (13,673 pairs, 7,545 distinct X):
reference only, not complete truth. Correspondence by sample, contig, interval (|delta| <= 10 bp = same interval;
reciprocal overlap >= 0.5 = boundary differs) and relative orientation; no family / pair IDs (`dna_validate.py`).

| class | v0.1.2 | v0.2 |
|---|---|---|
| recovered, same interval | 5,794 | 11,098 |
| one X recovered; pair deferred | 5,029 | 0 |
| one X recovered; pair boundary unresolved | 2,758 | 2,483 |
| not recovered: boundary unresolved | 81 | 81 |
| correspondence uncertain | 8 | 8 |
| no seed | 3 | 3 |
| protein strict X found in a DNA strict pass row | 7,371 / 7,545 | 7,371 / 7,545 |

Controls (protein-entry rows judged only_one_side_different between these assemblies): 17 / 17 not passed in the DNA
entry. DNA strict pass pairs without a protein-entry strict counterpart (kept; DNA-only candidates): 295,549.

v0.1.2 -> v0.2 concordance on the same occurrences (`concordance_v0.1.2_to_v0.2.json`): 0 of 7,211,086 source pairs
missing; 7,029,000 identical; 182,086 changed (wider grid windows / per-pair search space; all changed X re-checked);
3,282,617 only in v0.2 (formerly deferred). v0.2.0 vs v0.2.1 (`determinism_v0.2.0_vs_v0.2.1.json`): 10,493,701 of
10,493,703 identical; 2 changed between two unresolved states.

Scope: a segment is discovered only if it is a seed in at least one seed chunk or similar to a library sequence; X is
the shared extent within +/- 3 kb (grid-widened) windows, so pairs whose shared block runs past the window stay
unresolved. Not an exhaustive census, and not a claim of mobility.

## v0.2.3 / v0.2.4 (2026-11-14): review fixes, corrected pair counts, context-sample pairing

Code review (own review + Codex) fixes are in v0.2.3; the 50-assembly validation was rerun end to end (`v0.2.3/`).

**Correction of the pair counts above.** Before v0.2.3 a pair's ID included the seed order and the library-orientation
strands, so one pair of X intervals found through several library sequences was checked and counted several times.
With the ID taken from the two intervals and their relative orientation only:

| 700 <= X < 2500, both checks pass | v0.2.2 (duplicates) | v0.2.3 |
|---|---|---|
| pairs | 307,257 | **175,473** |
| distinct X intervals | 38,708 | 38,682 |
| same / cross assembly | 22,740 / 284,517 | 13,030 / 162,443 |
| pairs sent to the checks | 572,344 | 353,949 |
| cores with >= 2 / 4 / 6 / 8 contexts (any) | 525 / 48 / 34 / 25 | 524 / 48 / 34 / 25 |
| protein strict pairs recovered, same interval | 11,040 / 13,550 | 11,040 / 13,550 |

Intervals, context tiers and recovery are unchanged; only the pair count was inflated.

**Context-sample pairing (v0.2.4, `v0.2.4_context_sample/`).** At family scale all pairs cannot be formed
(Microbacteriaceae, 3,266 assemblies: 172 billion pairs over 394,613 library sequences). Each library sequence with
more than 16 occurrences gets a priority chunk: all pairs among 16 sampled occurrences (up to 8 copies from its top
assembly, then one per assembly in an assembly order shared by all library sequences). It is a pairing priority, not
a rule: all other pairs stay deferred and are counted as pending, and a core below a tier with pending pairs is
`insufficient_evidence`, not below the tier. On the validation set, against all pairs (v0.2.3):

| | all pairs | context sample |
|---|---|---|
| source pairs aligned | 10,493,703 | 770,449 |
| pairs sent to the checks | 353,949 | 44,990 |
| strict 700-2500 pairs | 175,473 | 3,924 |
| cores >= 4 / 6 / 8 contexts (any) | 48 / 34 / 25 | 43 / 24 / 11 |
| cores >= 4 / 6 / 8 (within one assembly) | 41 / 29 / 18 | 35 / 19 / 5 |
| cores >= 4 / 6 (distinct assemblies) | 16 / 7 | 9 / 3 |
| cores marked insufficient_evidence | 0 | 216 of 466 |

The sample finds 43 of the 48 main-list cores with 7% of the alignments and 13% of the checks; it undercounts the
higher tiers (a 16-occurrence sample cannot show every context of a high-copy core). Tier counts from a sampled run
are lower bounds.
