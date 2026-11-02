# Read-based site search: pilot design (draft, not started)

Rounds 1–5 took positions already in assemblies and checked them with reads. This stage asks a different
question: **given a known core, do the reads hold it at a position that the assembly catalog does not list?**

```
known core → find the core in reads → extract both flanks from each read → merge observations of one position
           → compare with the assembly position catalog → report candidate new sites
```

- **No empty site needed:** a position is reported from reads that carry the core with both flanks placed. It
  does not need an assembly without the core at that position.
- **One read hit is not one site:** hits are merged by molecule, then by position. Only the merged position is a
  candidate, and only with enough independent molecules.

## Pilot scope

- **Assemblies:** about 10, chosen for read quality, not for expected results.
  - Prefer molecule-level counting (ONT `molecules`, HiFi), long reads, and depth of at least 30 molecules.
  - Exclude runs whose reads are not raw sequencing output, e.g. SRR7743797 (submitted as `aligned_readsV3.bam`,
    reads in reference order and probably trimmed to aligned segments; see the round 5 notes).
- **Cores:** three groups, fixed before the search runs.
  - **Test cores:** cores with a supported A/B pair whose occurrences include a pilot assembly, of several product
    classes (transposase, recombinase, hypothetical, other).
  - **Negative control:** single-copy conserved genes (e.g. ribosomal proteins, `rpoB`). New positions found for
    these measure false positives: chimeric reads, mis-placed flanks.
  - **Positive control:** catalog positions already supported in rounds 1–5 in the same assemblies. Recovering
    them measures sensitivity.

## Steps

1. **Find the core in reads.**
   - Query: the core DNA of each occurrence (the tested element interval).
   - Align the reads of a pilot run to the core panel with minimap2, using the read-type preset.
   - Keep hits that cover most of the core. Record the read coordinates, strand and identity.
   - Hits that cover only part of the core are kept as a separate class, not dropped silently.
2. **Extract both flanks** from each hit: the read sequence beyond each end of the core hit.
   - A hit is *two-sided* when both flanks are at least L long (L to be fixed; ≥ 1 kb proposed, well above the
     caller's 300 bp anchor).
   - One-sided hits are recorded but do not by themselves make a position.
3. **Place each flank** on the sample's own assembly, then on the pool if it does not place there.
   - Keep every placement and its score.
   - A flank with two near-equal placements is `flank_ambiguous`. The same genome-wide competition rule as the
     current caller applies.
4. **Merge observations.**
   - First collapse reads of one molecule (same molecule grouping as the check stage).
   - Then group molecules whose U and D placements agree within a tolerance and have the same orientation.
   - The unit counted per position is molecules, or assumed molecules / read records where that is all a run
     allows.
5. **Compare with the catalog** (all occurrences of that core and its homologs in the pool), giving one class per
   merged position:
   - `catalog_position`: matches a listed occurrence in this assembly;
   - `catalog_position_other_assembly`: flanks place at a listed occurrence of another assembly only;
   - `unlisted_in_own_assembly`: both flanks place in this assembly, adjacent or close, at a position with no copy
     of the core. Reads carry the core at a position the assembly does not show;
   - `flanks_unplaced` / `flank_ambiguous` / `one_sided`: not a position.
6. **Report candidate new sites.** A candidate is an `unlisted_in_own_assembly` or
   `catalog_position_other_assembly` position with at least N independent molecules (N to be fixed; 3 like the
   current MIN_SPAN proposed). Each candidate lists:
   - molecules, read records, flank lengths and identities;
   - whether any molecule at the same place lacks the core (reported, not required);
   - the run and its count mode.

## What the pilot measures before any claim

- **Sensitivity:** recovery of the positive-control catalog positions.
- **False positives:** new positions for the negative-control cores, and for the test cores in assemblies where
  the core is absent.
- **Effect of thresholds:** how candidate counts change with L, N and the identity cutoff. The thresholds are
  fixed before the full run, not tuned to its results.

## How results will be worded

- **Same reads as the assembly:** the reads are mostly the ones the assembly was built from. An
  `unlisted_in_own_assembly` position is a position present in some molecules of the sequenced culture, or one the
  assembly collapsed. It is not an independent insertion event, and not mobility.
- **Not a single-read count:** a candidate is a merged position, never one read hit.
- **No result is not absence:** a core not found in reads means "not found here", not absent.

## Open decisions (for review before the pilot runs)

- Which ~10 assemblies and which cores.
- L (flank length), N (molecules), the identity and coverage cutoffs for core hits, and the placement tolerance.
- Whether the 104 assemblies left out of the read-check download plan should be considered for this stage. That
  plan scored assemblies only by known catalog positions. They may still matter for finding new sites.
