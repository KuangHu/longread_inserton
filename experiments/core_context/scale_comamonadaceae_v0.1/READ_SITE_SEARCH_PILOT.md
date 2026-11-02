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

## Pilot results (v0.1, 2026-11-01)

Files: `read_site_search_pilot/` (`summary.json` with the full grid, `positions.tsv.gz` with every merged position,
and one table per report group). Job ids in `read_site_search_pilot/jobs.txt`; panel summary in `panel.json`.

**What ran:** the 11 frozen samples, 12 runs, every read of every run.
- The last read with a panel hit lies within the final 0.1% of each read file, so the search reached the end of
  every file.
- The panel holds 1,083 cores: 18,797 distinct sequences from 35,666 occurrences, 17.6 Mb.

**Rules changed after the first run, before these results** (all found with the negative controls; no threshold
changed):
1. The two read strands of one position were first counted as two positions. They are now merged in assembly
   orientation (`canonical`).
2. A core hit that covers only part of a gene moves a flank anchor and made a second, offset position. Such a
   position overlapping a known one is now `shifted_known_position` and is not counted as a candidate.
3. "Does the assembly between the flanks hold the core" first used only the panel sequence, which can come from
   another genus. It now also uses the position's own reads (up to 3).
4. Positions where the assembly holds the read's core sequence but the catalog lists nothing are split by the
   annotated gene there (`annotate_omissions`). Only a gene of the core's own family, or no gene, counts as a
   catalog omission.
5. A junction joining the end and the start of one contig is `contig_end_join`.

### Main setting: flanks 1 kb, at least 3 independent molecules

| report group | class | positions |
|---|---|---:|
| known positions recovered | `known_catalog_position` | 533 |
| | `known_catalog_position_other_core` | 526 |
| | `shifted_known_position` (not candidates) | 93 |
| catalog omissions | `catalog_omission:same_family` | 0 |
| | `catalog_omission:no_annotated_cds` | 6 |
| homologs the catalog does not cover | `homolog_not_in_catalog:same_homology_group` | 295 |
| | `homolog_not_in_catalog:other_family` | 59 |
| structures the assembly does not show | `not_in_assembly:core_absent` | 4 |
| | `not_in_assembly:junction` | 0 |
| assembly breaks | `contig_end_join` | 5 |
| unplaceable contexts | `unplaceable_context:*` | 150 |

- **Recovery of known positions:**
  - Of the catalog positions of panel cores in these samples that read-check round 5 supported, 329 / 366 (90%)
    are recovered. Of the other catalog positions, 204 / 227 are recovered.
  - The 37 round-5-supported misses:
    - 3 were found with fewer than 3 molecules.
    - 5 have only unplaceable contexts nearby (flanks in repeats).
    - 1 has nothing placed near it.
    - 28 have a merged position overlapping or next to them that was matched to a different catalog occurrence.
      These are mostly two-gene occurrences and partial hits: a position is matched to an occurrence only if the
      occurrence lies wholly within it. Matching by overlap is a change for v0.2. It is not made here.
- **Negative control:**
  - Of the 55 (sample, gene) combinations searched, 47 give exactly one position and none gives two or more.
  - 8 give none. In those, the sample's own copy was not in the panel and only a distant homolog was.
  - The positions sit at the genes themselves: 30 as the control gene, and 15 under the catalog entry of a reviewed
    core that is the same gene.
  - The other 2 are the gene in a sample whose own copy was not in the control catalog. They are labelled
    `homolog_not_in_catalog:other_family` only because a control gene has no protein-family key to match.
- **Catalog omissions:** no position holds a gene of the core's own family that the catalog misses. The 6 left have
  no annotated gene there.
- **Homologs the catalog does not cover:** 295 + 59 positions hold the read's core sequence in a gene of another
  protein family, which the reviewed catalog does not cover. The 295 are in the core's homology group.
  - They are DNA-level homologs found by an 80%-coverage DNA hit, not catalog gaps of the core and not new sites.
  - Common products: electron transfer flavoprotein α, tripartite tricarboxylate transporters, catalase-peroxidase,
    dehydrogenases.
- **Assembly breaks:** the 5 `contig_end_join` positions are all in SAMEA115676957, a contig-level assembly. The reads
  join the end of NZ_CAXUQU020000001.1 (6,888,869 bp) or NZ_CAXUQU020000002.1 (487,881 bp) to the same contig's
  start, through the core. The assembly breaks inside the core there. This is compatible with closed replicons; the
  topology is not established here.
- **Candidate structures not shown by the assembly** (`not_in_assembly:core_absent`, all 4):

| sample | core | contig : site | molecules with the core | molecules spanning the empty site | flank overlap |
|---|---|---|---:|---:|---:|
| SAMN30525180 (ONT) | PF_a5c0db98e215 | NZ_JAOCIM010000001.1 : 3,553,051–3,553,059 | 41 | 16 | 8 bp |
| SAMN11265848 (ONT) | PF_104f0ed0dd1b | NZ_CM195847.1 : 3,064,452–3,064,457 | 24 | 177 | 5 bp |
| SAMN15738716 (HiFi) | PF_7d1bf17c86e7 | NZ_JARXUY010000001.1 : 3,158,294–3,158,304 | 10 | 53 | 10 bp |
| SAMN30525180 (ONT) | PF_a5c0db98e215 | NZ_JAOCIM010000001.1 : 1,653,143–1,653,151 | 9 | 48 | 8 bp |

  - **What the reads show:** at each site the assembly shows the empty site, and some molecules of the sequenced
    DNA carry the core there with unique 1 kb flanks on both sides. Other molecules span the empty site. Both forms
    are present.
  - **Flank overlap:** the two flank anchors overlap by 5–10 bp, which fits a short duplication of the target
    sequence at the junctions. That is an observation, not a test.
  - **The cores:** each is multi-copy in its sample (11–23 recovered known positions), and each class count above
    is unchanged from 2 to 5 molecules and from 500 bp to 2 kb flanks.
  - **Wording:** these are candidate structures present in part of the sequenced DNA. They are not independent
    insertion events, and mobility is not shown.
- **Unplaceable contexts:** 150 merged contexts, by side:
  - 91 have one flank placed and the other ambiguous;
  - 40 have one flank placed and the other unplaced;
  - 19 have both flanks unplaced or ambiguous.

  They are not called new and not called absent.

### Threshold sensitivity (core panel; recovery = catalog positions recovered)

| flank | N | known recovered (round-5 supported) | known recovered (other) | known (same / other core) | catalog omission | homolog not in catalog (group / other) | contig end join | core absent | junction | unplaceable | neg: genes with 1 / >1 position |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 500 | 1 | 326/366 | 201/227 | 600 / 717 | 7 | 435 / 78 | 7 | 15 | 56 | 1381 | 41 / 8 |
| 500 | 2 | 323/366 | 201/227 | 527 / 564 | 7 | 342 / 66 | 5 | 4 | 2 | 369 | 45 / 2 |
| 500 | 3 | 323/366 | 200/227 | 523 / 530 | 6 | 307 / 61 | 4 | 4 | 0 | 196 | 47 / 0 |
| 500 | 5 | 322/366 | 198/227 | 518 / 512 | 6 | 273 / 56 | 3 | 4 | 0 | 126 | 45 / 0 |
| 1000 | 1 | 332/366 | 206/227 | 657 / 746 | 7 | 445 / 84 | 11 | 18 | 69 | 1171 | 38 / 10 |
| 1000 | 2 | 329/366 | 206/227 | 542 / 566 | 6 | 332 / 63 | 6 | 4 | 3 | 296 | 45 / 2 |
| **1000** | **3** | **329/366** | **204/227** | **533 / 526** | **6** | **295 / 59** | **5** | **4** | **0** | **150** | **47 / 0** |
| 1000 | 5 | 329/366 | 203/227 | 531 / 508 | 6 | 269 / 52 | 4 | 4 | 0 | 99 | 44 / 0 |
| 2000 | 1 | 342/366 | 210/227 | 682 / 745 | 7 | 424 / 85 | 13 | 15 | 62 | 1160 | 34 / 9 |
| 2000 | 2 | 338/366 | 210/227 | 555 / 563 | 6 | 310 / 63 | 8 | 4 | 4 | 349 | 39 / 3 |
| 2000 | 3 | 338/366 | 207/227 | 547 / 524 | 6 | 281 / 58 | 7 | 4 | 0 | 134 | 42 / 0 |
| 2000 | 5 | 338/366 | 205/227 | 543 / 506 | 6 | 254 / 54 | 7 | 4 | 0 | 52 | 40 / 0 |

- **One molecule is not enough:** at N = 1, junctions (56–69) and core-absent positions (15–18) appear, and
  negative-control genes get extra positions (8–10). All of these vanish by N = 3. This is the expected
  signature of chimeric or mis-placed single reads, and it is why one read hit is never counted as a site.
- **The candidates are stable:** the 4 core-absent positions are the same at every flank length from N = 2 up.
- **Flank length trades placement for reach:** longer flanks recover a few more known positions (329 → 338 at
  N = 3) but leave fewer molecules with both flanks; shorter flanks leave more contexts unplaceable.

### Not done in the pilot

- Placing unplaced or ambiguous flanks on other assemblies of the pool.
- Matching positions to catalog occurrences by overlap rather than containment (the 28 recovery misses).
- Any assembly outside the frozen 11; no downloads.
