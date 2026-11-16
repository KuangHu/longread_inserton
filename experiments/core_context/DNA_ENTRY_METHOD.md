# DNA-first search for one core in several different contexts (method, v0.2.5)

## What it looks for

A DNA segment X ("core") that occurs at two or more positions where the 1,000 bp before it (U) and the 1,000 bp
after it (D) are both different:

    position A:  U_A - X  - D_A
    position B:  U_B - X' - D_B        X and X' similar, U_A / U_B different, D_A / D_B different

The result is a structural pattern in assemblies. It is not a call of mobility, of an insertion event or of an
empty site. No protein, gene annotation, product name, gap caller, site grouping or read evidence is used to admit
a candidate; read support and annotation are separate fields added later.

## Input

Assembly FASTA files listed in a manifest with their SHA-256. Every file is hashed before it is read; a mismatch is
recorded and the assembly is not used. Analysis jobs never download.

## Steps

1. **Seeds.** Assemblies of one family are split into chunks of 75 (related genomes together). Within a chunk,
   all assemblies are aligned against each other (minimap2 `asm20`, all chains kept). An alignment of 100-10,000 bp
   on both sides with identity >= 80% is a seed: a repeated or shared segment. Seed ends are where the aligner
   stopped, not boundaries. Seeds are clustered (90% identity, 80% coverage, containment) into representatives.
2. **Library.** Representatives of all chunks are clustered again into one library of sequences for the family.
3. **Occurrences.** Every library sequence is searched in every assembly (minimap2 `asm20`, no copy cap; hits
   >= 100 bp and >= 80% identity). Overlapping hits of one library sequence on one contig are one occurrence.
4. **Pairing.** Occurrences of the same library sequence are paired, within one assembly and across assemblies.
   All pairs cannot be formed at family scale (Microbacteriaceae: 172 billion), so each library sequence with more
   than 16 occurrences gets a context sample: all pairs among 16 occurrences (up to 8 copies from the assembly with
   most copies, then one per assembly in an order shared by all library sequences). The sample is a priority, not
   a rule: every other pair stays deferred and is counted.
5. **Pair-specific X.** For each pair, the two windows (occurrence +/- 3,000 bp) are aligned with blastn and the
   shared block around the seed is chained and extended. X is the seed plus the shared extension on each side, and
   it is resolved only if both ends lie inside the windows and contigs in both copies. X is never trimmed.
   Identical alignments reached through different library sequences are computed once.
6. **Target length.** Pairs with 700 <= X < 2,500 bp on both sides go to the checks; others are kept as
   out of target. Duplicate pairs (same two intervals) are kept once.
7. **Two checks** (frozen thresholds, one blastn per pair):
   - Check 1, source: U, X and D are extracted again from the source assembly, whose hash must match.
   - Check 2, structure: X / X' similar (identity >= 90%, coverage >= 90% of both); U_A / U_B and D_A / D_B both
     different (<= 20% of each flank covered by >= 80%-identity alignments). Flanks with fewer than 800 usable bases,
     flanks that are similar or uncertain, swapped or rearranged contexts and pairs that reuse one physical
     stretch do not pass. Missing context is never counted as different.
8. **Distinct contexts per core.** X intervals overlapping >= 50% on one contig are one position. Positions linked
   by a passing pair form a core. The number of contexts of a core is the largest set of its positions in which
   every pair passes (so three copies that share one context count once). Counted three ways: any positions,
   all within one assembly, all in different assemblies.

## Tiers (project rule)

2-3 contexts exploratory; >= 4 main list; >= 6 priority; >= 8 showcase. A core below a tier whose library sequences
still have deferred pairs is `insufficient_evidence`, not below the tier.

## What the numbers mean and do not mean

- Counts from context-sample runs are lower bounds: most pairs are deferred, and a 16-occurrence sample cannot show
  every context of a high-copy core. On the 50-assembly validation set the sample found 43 cores with >= 4 contexts
  against 48 with all pairs, and 11 against 25 with >= 8.
- A segment is found only if two copies meet in one seed chunk or it matches a library sequence. The search is not
  exhaustive.
- At family scale a conserved gene whose flanks have diverged below 80% identity passes both checks across
  assemblies. The within-one-assembly count is the more specific one; the across-assembly count includes such cases.
- Families are run separately; cores are not compared between families.

## Running it

- `dna_fetch.sbatch <family> <taxon> <root> [queue.tsv]`: the only job that downloads; one CPU on the shared
  partition `lr7` (charged one CPU), at most 24 h; writes the verified manifest and appends the family to the queue.
- `dna_queue.sh run <root> <queue.tsv> [parallel]`: per family, in its own work directory: input check -> run
  (`dna_run.sh`) -> verify -> archive with checksums -> delete the listed re-derivable intermediates, recording the
  space freed. `dna_queue.sh status <root> <queue.tsv>` prints each family's stage and tier counts.
- `lr6` charges a whole node per job, so heavy steps use all 32 CPUs and small tasks share a node job.
- Results per family: `<root>/archive/<family>/` (report, tier tables `core_contexts/`, strict pair table) and
  `<root>/<family>/work/final_<family>/` (all rows).

## Results so far (context sample; lower bounds)

| family | assemblies | cores >= 2 | >= 4 / 6 / 8 (any) | >= 4 / 6 / 8 (within one assembly) | >= 4 / 6 / 8 (different assemblies) |
|---|---|---|---|---|---|
| Microbacteriaceae | 3,266 | 10,459 | 508 / 159 / 70 | 179 / 71 / 28 | 341 / 63 / 15 |
| Comamonadaceae | 1,686 | 11,765 | 701 / 259 / 136 | 317 / 137 / 42 | 428 / 83 / 20 |
| Micrococcaceae | 2,972 | 12,300 | 922 / 310 / 138 | 322 / 112 / 35 | 610 / 139 / 41 |
| Lachnospiraceae | 2,970 | 11,162 | 889 / 385 / 208 | 453 / 228 / 78 | 536 / 123 / 33 |
| Actinomycetaceae | 1,133 | 1,686 | 137 / 58 / 21 | 91 / 33 / 9 | 56 / 10 / 1 |
| Metamycoplasmataceae | 1,583 | 603 | 79 / 51 / 37 | 66 / 39 / 14 | 38 / 18 / 9 |
| Eggerthellaceae | 315 | 781 | 28 / 6 / 4 | 16 / 4 / 2 | 16 / 1 / 1 |

Flavobacteriaceae, Nitrobacteraceae and Erwiniaceae were still running on 2026-11-15.
