# Read-based site search: pilot design (v0.1)

Rounds 1–5 took positions already in assemblies and checked them with reads. This stage asks a different
question: **given a known core, do the reads hold it at a position that the assembly or its catalog does not
show?**

```
known core → find homologous cores in all reads → take both flanks from each read → merge observations of one
position by molecule → compare with the assembly and its position catalog → report by class
```

- **No empty site is needed.** A position is reported from reads that carry the core with both flanks placed.
- **One read hit is not one site.** Hits are merged by molecule, then by position; only merged positions are
  counted.
- **All qualifying reads are searched.** Every read of every verified run of a pilot sample is aligned to the core
  panel, not only reads near catalog positions.
- **No new downloads.** The pilot uses reads already verified on disk.
- Code: `../read_site_search.py` (stages `freeze`, `panel`, `search`, `combine`) and `../read_site_search.sbatch`.

## Frozen pilot samples (`read_site_search_pilot/`, frozen before any search)

- **Rule** (`freeze.json`; every assembly considered is in `candidates_considered.tsv` with the reason it was not
  eligible):
  - **Eligible:** every checked run counted in molecules; read type ONT or PacBio HiFi/CCS; no run submitted as an
    aligned BAM; depth at least 30× (verified bases / assembly length).
  - **Per genus:** the eligible sample with the most candidate cores supported in read-check round 5 (ties: depth,
    then id). Genera whose best sample has at least 15 such cores.
  - **HiFi minimum:** then HiFi samples by the same score until at least 2 HiFi samples are included.
- **Result:** 11 samples from 10 genera, 9 ONT and 2 HiFi (`samples.tsv`, with run ids and SHA256 of the read
  files).

## Core panel and catalog

- **Test cores:** the 1,083 candidate cores with an A/B pair supported on both sides in round 5. Their panel holds
  the distinct DNA sequences of the coding span of every occurrence in the pool, in core orientation.
- **Catalog:** every reviewed occurrence (any core) in the pilot samples. A position is compared with the catalog
  of its own core and with that of any other core.

## Controls

- **Known positions (positive control):** catalog occurrences of panel cores in the pilot samples, split into those
  supported in round 5 and the rest. Recovery is the fraction found as `known_catalog_position` with enough
  molecules. It measures sensitivity, and its misses show where the search loses positions.
- **Single-copy genes (negative control):** per pilot sample, genes from a fixed product list (RNA polymerase β and
  β′, gyrase B, RecA, DNA polymerase III β, ribosomal proteins L2 and S3, Ffh). Each must occur exactly once in the
  annotation, complete and not pseudo, and must overlap no catalog occurrence. They go through the same search.
  - **Panel in practice (`panel/panel.json`):** 30 genes. These are RNA polymerase β and β′ in all 11 samples, L2
    in 3, S3 in 4 and DNA polymerase III β in 1.
  - **Products that gave no gene:** gyrase B, RecA and Ffh. Their annotated names differ from the list, or the
    gene overlaps a catalog occurrence. The list was not changed after the panel was built.
  - Expected: exactly one position each, `known_catalog_position`.
  - Any other position for them measures false positives from chimeric reads, mis-placed flanks or the merging
    rules.
- **Thresholds:** fixed before the run. The main setting is flanks of 1 kb and at least 3 independent molecules.
  Sensitivity is reported for flanks of 500 / 1,000 / 2,000 bp × 1 / 2 / 3 / 5 molecules. No threshold is
  tuned on the results.

## Steps

1. **Core hits:** minimap2 (read preset) of all reads against the panel.
   - A hit covers at least 80% of a panel sequence.
   - Per read and core, overlapping hits keep the best one.
2. **Flanks:** exactly L bp of read sequence next to each end of the hit. A hit without L bp on both sides is
   `one_sided` at that L.
3. **Placement:** each flank is aligned to the sample's own assembly.
   - `placed`: the best alignment covers at least 80% of the flank and no other alignment within 20 AS
     (MIN_DELTA_AS) places it elsewhere (beyond 50 bp, another contig or strand).
   - `ambiguous`: such a competing placement exists.
   - `unplaced`: no alignment covers enough of the flank.
4. **Merge:**
   - Reads of one molecule count once (the same grouping as the read check).
   - Molecules whose two flank anchors agree within 50 bp form one position.
5. **Classes** (see the docstring of `read_site_search.py` for the exact rules):

| report group | class | meaning |
|---|---|---|
| known positions recovered | `known_catalog_position` | a catalog occurrence of this core lies between the two flanks |
| | `known_catalog_position_other_core` | a catalog occurrence of another (homologous) core lies there |
| catalog omissions | `catalog_omission` | the assembly sequence between the flanks contains the core, but no catalog occurrence is listed there |
| structures the assembly does not show | `not_in_assembly:core_absent` | the flanks are adjacent in the assembly; the reads carry the core at a junction the assembly shows without it |
| | `not_in_assembly:other_sequence` | sequence of core-like length between the flanks that does not contain the core |
| | `not_in_assembly:junction` | flanks on different contigs, opposite strands or far apart |
| unplaceable contexts | `unplaceable_context:<left>\|<right>` | at least one flank unplaced or ambiguous, reported by side |

## How results are worded

- **Same reads as the assembly:** the reads are mostly the ones the assembly was built from.
  - A `not_in_assembly` position means some molecules of the sequenced culture carry the core at a junction the
    assembly does not show. It can be a minority structure in the culture, a structure the assembler collapsed or
    resolved differently, or chimeric molecules.
  - It is not an independent insertion event, not mobility, and not a new site in any population sense. It is
    called a **candidate structure not shown by the assembly**.
- **Catalog omissions** are gaps in the catalog (annotation or review filters), not new sites.
- **Unplaceable contexts** are flank sequences the own assembly does not hold uniquely. They are not called absent
  from the assembly, and they are not called new.
- **Merged positions only:** every count is of merged positions with their molecules, never of read hits.
- **Not found is not absent:** a known position that is not recovered is a sensitivity miss to explain, not
  evidence that the structure is absent.

## Not in the pilot

- Placing unplaced flanks on other assemblies of the pool.
- Assemblies outside the frozen 11, including the 104 left out of the read-check downloads.
