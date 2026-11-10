# U–insertion–D comparison v0.1 (`../../ucd_compare.py`), small-scale run

**Question:** for two occurrences `U_A – X – D_A` and `U_B – X′ – D_B`, is X similar to X′ while U and D both differ?

- The answer is a structural pattern only. It says nothing about mobility, and no empty site is required.
- No gap caller, S6 site group or read-check result is used to admit a pair. Their labels are carried as columns for
  comparison.
- Nothing is re-annotated, downloaded or searched in reads. Sequences come straight from the manifest assemblies.

## Run

```bash
D=scale_comamonadaceae_v0.1
python ucd_compare.py prepare --out $D/ucd_compare_v0.1 \
  --samples $D/read_site_search_pilot/samples.tsv,$D/read_site_search_validation1/samples.tsv,$D/read_site_search_validation2/samples.tsv
python ucd_compare.py run --out $D/ucd_compare_v0.1 --work /global/scratch/users/kh36969/longread_callout/core_context/ucd_compare_v0.1 --threads 16
python ucd_compare.py figures --out $D/ucd_compare_v0.1 --work /global/scratch/users/kh36969/longread_callout/core_context/ucd_compare_v0.1
```

SLURM job 27590448 ran on one lr6 node with 16 CPUs: 4 min wall, 102 MB RSS (`jobs.txt`).

## Rules (v0.1 working definition; `CONFIG` in the code, copied into `summary.json`)

Rules were set before the run and not tuned on its results.

| step | rule |
|---|---|
| sequences | X = interval; U / D = 1,000 bp before / after X, all in the insertion's orientation (a minus-strand occurrence is reverse-complemented) |
| insertion similar | reliable HSPs (≥ 90% identity, ≥ 30 bp, e ≤ 1e-5) cover ≥ 90% of X and of X′ on one strand; some reliable HSP but less coverage → `insertion_partial`; none → `insertion_dissimilar` |
| orientation | if X′ matches X on the minus strand, B is flipped (U_B′ = rc D_B, D_B′ = rc U_B) before contexts are compared |
| context comparison | non-redundant union of reliable HSPs (≥ 80%, ≥ 30 bp) on each flank: **similar** both ≥ 80%; **different** both ≤ 20%; **uncertain** otherwise |
| checked comparisons | U/U and D/D same strand; U/U and D/D inverted; U_A/D_B and D_A/U_B on both strands |
| insufficient context | a flank with < 800 bp that are not N and not in low-complexity windows (64 bp, 3-mer entropy < 3 bits); contig edges shorten flanks |
| shared physical flank | same assembly and contig, and a flank of one occurrence overlaps a flank or the insertion of the other |

**Verdict:** the first rule that applies wins, and every applicable reason is listed in `reasons`.
1. `insertion_dissimilar` / `insertion_partial`
2. `shared_physical_context`
3. `insufficient_context`
4. `swapped_or_rearranged_context`: an inverted or cross-side comparison is similar
5. `same_context`
6. `only_one_side_different`
7. `ucd_supported`: U different, D different, and every inverted or cross-side comparison different
8. `uncertain_context`: everything else

**Interval kinds** are always kept apart:
- `cds_core` is the catalog CDS of the core.
- `candidate_insertion` is the coding span plus the pair's own shared extension, where both boundaries were resolved
  by the boundary review. It is pair-specific.
- A CDS is never reported as a complete insertion.

**Invariance:** pair order is canonical (A = the smaller occurrence ID). Tests check that swapping A and B,
reverse-complementing either side, or shuffling the input does not change the verdict.

## Input (`input_pairs.tsv.gz`, `prepare.json`)

- **918 pairs:**
  - 901 read-check round-5 pairs whose two samples are both among the 24 development / validation samples
    (pilot 11, validation 1–2 13): tier A (same assembly) 563, tier B (cross assembly) 51, and 287 "repeated
    context" controls of the earlier step;
  - 17 named regression cases (UCD01–12, UCDR01, 03–06).
- UCDR02 is absent because at pool scale its copies are two core families, so there is no core pair.
- Labels copied for viewing only: tier, product, control class, core match, boundary status, estimated span,
  read-check position / pair class, and the earlier S6 context class (`reference_context_class`).

## Results (`ucd_pairs.tsv.gz`, one row per pair and interval kind; `summary.json`)

Counts are pairs; occurrences, samples and cores are counted separately. Pairs are not independent events.

| verdict | cds_core pairs / occ. / samples / cores | candidate_insertion pairs / occ. / samples / cores |
|---|---|---|
| ucd_supported | 442 / 419 / 32 / 120 | 573 / 515 / 33 / 157 |
| uncertain_context | 134 / 133 / 24 / 55 | 9 / 15 / 5 / 7 |
| only_one_side_different | 51 / 66 / 11 / 28 | 0 |
| same_context | 190 / 366 / 15 / 173 | 0 |
| swapped_or_rearranged_context | 1 / 2 / 1 / 1 | 1 / 2 / 1 / 1 |
| shared_physical_context | 13 / 22 / 9 / 9 | 0 |
| insufficient_context | 1 / 2 / 1 / 1 | 1 / 2 / 1 / 1 |
| insertion_partial | 3 / 6 / 2 / 3 | 2 / 4 / 1 / 2 |
| insertion_dissimilar | 83 / 139 / 21 / 56 | 44 / 74 / 20 / 32 |
| interval_not_available | — | 288 / 535 / 23 / 234 (287 controls have no candidate interval; UCDR06 has an undetermined extension) |

"Samples" counts both sides' samples of cross-assembly pairs.

**By source:**
- **Named cases:** all 17 are `ucd_supported` on `cds_core` or `candidate_insertion` (UCDR06 has no candidate
  interval).
  - On `cds_core`, 3 are `uncertain_context` (UCD08, 09, 12): the shared sequence continues past the CDS into the
    flank.
  - On the candidate interval, which ends where the shared sequence ends, all three are `ucd_supported`.
- **Tier A** (same assembly, 563 pairs):
  - `candidate_insertion`: 503 ucd_supported, 8 uncertain, 25 insertion_dissimilar, 2 partial, 1 rearranged,
    1 insufficient (+23 conserved-paralog controls, all ucd_supported).
  - `cds_core`: 378 ucd_supported, 105 uncertain, 35 one side different. The CDS boundary often lies inside the
    shared extension, so its flanks start with shared sequence.
- **Tier B** (cross assembly, 51 pairs): `candidate_insertion` 31 ucd_supported, 19 insertion_dissimilar,
  1 uncertain.
- **Earlier "repeated context" controls** (287, `cds_core` only):
  - 190 same_context, 13 shared physical context, 48 insertion_dissimilar (a protein family whose copies differ at
    DNA level), 3 partial, 11 one side different, 16 uncertain;
  - **6 ucd_supported:** copies 3.6–10 kb apart on one contig. Shared DNA ends 50–180 bp outside the CDS (5–18% of
    each 1 kb flank, at ~100% identity), and the rest of each flank differs.
  - The earlier step judged "repeated context" from neighbouring host genes at protein level. At DNA level these
    flanks are different, so the two definitions disagree here.
  - These pairs are listed as such, not removed.
- **Conserved paralogs** (EF-Tu and similar, `control_class=conserved_paralog_control`): 23 tier-A pairs are
  `ucd_supported` and are kept with that label, as the rule says.

**Insertion dissimilar despite one core family:** these are protein-level families (`homolog_full_coverage`, 35 on
`cds_core`) whose two copies are < 90% DNA identity, or whose reliable HSPs cover < 90%.

## Figures (`figures/`, SVG; `figures/index.json`)

One example per requested pattern, chosen by rule (a named case first, then the highest insertion identity):

| figure | pair | what it shows |
|---|---|---|
| clear pass | UCD02 (`candidate_insertion`) | X/X′ 100% over 824 bp; no reliable HSP in any of the 8 context comparisons |
| same context | VC_00a85a0dba74 (`cds_core`) | both flanks similar on the same side |
| only one side different | VP_05b4be5ebd7f (`cds_core`) | one side similar, the other different |
| orientation / swapped | VP_e90421fa9b5b (`candidate_insertion`) | U/U and D/D share nothing; U_A~rc(D_B) and D_A~rc(U_B) cover 100% both ways: X inverted in place |
| insufficient context | VP_fa59b4a0a82d (`candidate_insertion`) | a flank with < 800 usable bp |

Each figure shows the measured identity and coverage per comparison. Blocks under a flank are the reliable HSPs:
dark for same side and same strand, red for inverted or cross side.

## Evidence

Every comparison's oriented sequences (U, X, D of both sides) and all blastn HSPs are written to scratch. Paths, sizes
and sha256 are in `summary.json`:
- `evidence_sequences.fa.gz` (1.2 MB);
- `evidence_hsps.tsv.gz` (57 KB).

## Scale-up estimate

- **Cost:** 4.1 CPU-s per pair (both interval kinds, 4 blastn calls), measured on lr6 (238 s wall × 16 CPUs /
  918 pairs). Memory is negligible (102 MB). Evidence is about 1.3 KB per pair.

| set | pairs | CPU-h | wall on lr6 | evidence |
|---|---|---|---|---|
| this run | 918 | 1 | 4 min, 1 node | 1.2 MB |
| `match_both_lt2500` (both elements < 2.5 kb) | 95,328 | ≈ 110 | ≈ 7 h on 1 node of 16 CPUs, or ≈ 1 h as 8 array tasks | ≈ 125 MB |
| every core pair of the boundary review | 3,894,917 | ≈ 4,400 | ≈ 280 node-h | ≈ 5 GB |

- **Bottleneck:** one blastn process per pair. Batching many pairs per blastn call would cut the cost for the full
  set; that is not part of v0.1.
- `prepare` would need the pairs of the chosen set in the same columns (coordinates come from the boundary-review
  pair shards and the occurrence shards, as for the named cases here).
