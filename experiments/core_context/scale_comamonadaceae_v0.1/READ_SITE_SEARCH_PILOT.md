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
      **Corrected in v0.2 (below):** with matching by core relation, orientation and coverage, recovery is
      unchanged (329 / 366). These misses are not a matching problem; see the v0.2 recovery misses.
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

## Pilot v0.2 (2026-11-01): rule fixes, same 11 samples, same thresholds

Requested after review of v0.1. Panel, samples, panel PAFs and all thresholds are unchanged. Only the per-sample
search after the panel alignment and the combine step were rerun. Jobs are in `read_site_search_pilot/jobs.txt`.
The v0.1 per-sample tables are kept in the work directory under `v0.1_backup/`; the v0.1 combined tables are in git
at 8af82d6. Tests in `../test_read_site_search.py` include one counterexample per fix.

### What changed

1. **Flank inner end (`placement`):**
   - v0.1 took the inner coordinate from the end of the best flank alignment even when that alignment stopped short
     of the core.
   - A flank is now `placed` only if its best alignment reaches the inner end within `INNER_CLIP` = 20 read
     bases. This value was fixed before the rerun.
   - Otherwise the flank is `inner_end_unaligned`, and the context is reported as `unplaceable_context:*`.
   - Of the 217,223 aligned 1 kb flanks:

     | bases short of the inner end | share of flanks |
     |---|---:|
     | 0 | 84.9% |
     | 1–20 | 13.2% |
     | 21–50 | 1.1% |
     | more than 50 | 0.8% |
2. **Catalog matching (`match_catalog`):** an occurrence matches only if all three conditions hold.
   - **Core relation:** the occurrence is of the same core, or of a homologous core that shares a protein family or
     homology group.
   - **Orientation:** the occurrence's strand agrees with the core's orientation in the assembly as the reads show
     it. For a homologous core, the shared family's sign converts between the two.
   - **Coverage:** the overlap is at least 80% of the occurrence and at least 80% of the interval.
   - An unrelated core's occurrence is only noted (`catalog_unmatched`), whatever its overlap.
   - A same or homologous occurrence that fails orientation or coverage gives `catalog_mismatch:<reason>`. That
     class is neither known nor a candidate.
   - Control genes are related to cores through their own protein families.
3. **`shifted_known_position` (`mark_shifted`):**
   - A position counts as a shifted copy of a known position only if the known position's core is the same or
     homologous, in a consistent orientation.
   - Two known positions of one core matched to the same occurrence are one position. The one with fewer molecules
     is marked shifted. The negative control found this case: a partial hit of a control gene matched the same
     occurrence as the full hit.
   - `shifted_from` names the known position.
4. **Empty-site molecules (`empty_site_reads` / `empty_site_counts`):** a record counts if two conditions hold.
   - One of its whole-assembly alignments spans the site ± 1 kb.
   - Its local CIGAR over the site ± 300 bp passes the read check's indel rule (`cigar_ok`, limit half the insert
     length).

   v0.1 counted read names and used only the whole alignment's length difference. Records are now merged into
   molecules. A molecule that carries the core and also spans the empty site is counted apart (`both_forms`).

### Main setting (flanks 1 kb, at least 3 molecules), v0.1 → v0.2

| report group | class | v0.1 | v0.2 |
|---|---|---:|---:|
| known positions | `known_catalog_position` | 533 | 530 |
| | `known_catalog_position_other_core` | 526 | 600 |
| | `catalog_mismatch:coverage` (new; not known, not a candidate) | – | 64 |
| | `shifted_known_position` (not counted) | 93 | 92 |
| catalog omissions | `catalog_omission:same_family` | 0 | 0 |
| | `catalog_omission:no_annotated_cds` | 6 | 6 |
| homologs not in the catalog | `homolog_not_in_catalog:same_homology_group` | 295 | 138 |
| | `homolog_not_in_catalog:other_family` | 59 | 73 |
| structures the assembly does not show | `not_in_assembly:core_absent` | 4 | **0** |
| | `not_in_assembly:junction` / `other_sequence` | 0 / 0 | 0 / 0 |
| assembly breaks | `contig_end_join` | 5 | 3 |
| unplaceable contexts | `unplaceable_context:*` | 150 | 390 |
| | of which with a flank `inner_end_unaligned` | – | 243 |

- **Recovery:** round-5-supported 329 / 366 (unchanged); others 201 / 227 (v0.1: 204).
- **Negative control:** 47 genes with exactly one position, 8 with none, none with two or more (as in v0.1).
- **Where the v0.1 positions went** (`read_site_search_pilot/v01_to_v02.tsv`; v0.1 "catalog omission" is
  counted before the family split):

  | v0.1 class | v0.2 class | positions |
  |---|---|---:|
  | catalog omission (360) | homolog not in the catalog | 185 |
  | | known position of a homologous core | 123 |
  | | shifted | 24 |
  | | `catalog_mismatch:coverage` | 21 |
  | | catalog omission | 6 |
  | known position of another core (526) | the same class | 446 |
  | | `catalog_mismatch:coverage` | 43 |
  | | shifted | 22 |
  | | homolog not in the catalog (the core in that occurrence is unrelated) | 15 |
  | known position of the same core (533) | the same class | 531 |

  - The rise in homologous-core matches is because a homolog's occurrence now matches with 80% reciprocal coverage
    instead of lying wholly inside within 50 bp.
  - The 15 positions that became "homolog not in the catalog" had matched an occurrence of an unrelated core.
- **`catalog_mismatch:coverage`:** all 64 overlap a homologous core's occurrence that covers less than 80% of the
  interval or of itself. This is typically a one-gene occurrence inside a two-gene core's interval, or the reverse.
  No same-core occurrence failed on orientation.
- **The 92 shifted positions:**
  - 77 are a second position of one core on the same homologous occurrence (rule 3).
  - 15 overlap a known position of a related core.
- **Sensitivity:** junctions and core-absent positions at N = 1 (68 junctions, 0 core-absent at 1 kb) vanish by
  N = 2–3 as before. The negative control has 2–3 genes with two positions at N = 1 and none from N = 3 at every
  flank length.

### The 4 v0.1 candidates

All 4 left `not_in_assembly:core_absent` because of rule 1.
- **Why:** at each site the flank alignments stop a consistent number of bases short of the core hit on both
  sides. The reads therefore carry more inserted sequence than the core hit: the core plus about 60–105 bp of end
  sequence that the assembly lacks.
- **v0.1's coordinates:** v0.1 placed the junctions where the assembly alignment stops. Under the strict rule,
  those coordinates are not the inner ends of the flanks.
- **The site is still known:** the search no longer places these positions, but the junction where the reads
  leave the assembly is still observed.
- **Strict test:** each site was therefore tested with the read check's own rules (`../candidate_evidence.py`;
  s3 v0.3 rules, thresholds unchanged).
  - **Alleles:**
    - P = U + the read's whole insert (core hit plus the unaligned ends, from the best-identity record) + D;
    - A = the assembly as it is.
  - **Reads:** every record with a whole-assembly alignment within 2 kb, plus the core-carrying records.
  - **Counting:** molecules, merged across runs.
- **Outputs** (`read_site_search_pilot/candidates/`):
  - `evidence_summary.tsv`;
  - `evidence_molecules.tsv`: one row per molecule, with its records' calls and nothing dropped;
  - one SVG per candidate, one row per molecule on U | insert | D.

| sample : site | v0.1 molecules (core / empty) | v0.2 search | unaligned ends L / R (median, range) | insert / core hit | strict P matched | strict A matched | conflict | unknown | empty-site molecules, v0.2 rule |
|---|---|---|---|---|---:|---:|---:|---:|---:|
| SAMN30525180 NZ_JAOCIM010000001.1:3,553,051–059 | 41 / 16 | 40 molecules, both flanks `inner_end_unaligned` | 69 (58–91) / 36 (30–54) | 1,334 / 1,228 | 46 | 25 | 0 | 58 | 22 (+2 indel) |
| SAMN30525180 NZ_JAOCIM010000001.1:1,653,143–151 | 9 / 48 | 9, both `inner_end_unaligned` | 69 (65–75) / 35 (33–39) | 1,328 / 1,226 | 10 | 90 | 0 | 47 | 78 |
| SAMN11265848 NZ_CM195847.1:3,064,452–457 | 24 / 177 | 23 both, 1 one side `inner_end_unaligned` | 21 (17–60) / 76 (68–94) | 1,176 / 1,079 | 30 | 302 | 0 | 1,028 | 197 |
| SAMN15738716 NZ_JARXUY010000001.1:3,158,294–304 (HiFi) | 10 / 53 | 10, one side `inner_end_unaligned` | 57 (56–58) / 7 (7–14) | 1,402 / 1,338 | 12 | 61 | 0 | 38 | 53 |

- **Which side is which:** L / R are in assembly orientation. The two SAMN30525180 sites show the same end lengths
  (69 / 35–36 bp) with the same core, which fits the same inserted unit at both.
- **Strict P / A counts:** the caller's own counts, per run, agree with the per-molecule table
  (`per_run_caller_counts`).
- **The read check's verdict, taking the assembly's allele (A) as the expected one:**
  - 3,553,051: `mixed` (46 P / 25 A);
  - SAMN15738716: `mixed` (12 / 61);
  - 1,653,143: `consistent` with A, with P at 10% of 100;
  - SAMN11265848: `consistent` with A, with P at 9% of 332.

  "Consistent" means only that P is at or under the read check's 10% rule. It does not mean P is absent: 10 and
  30 molecules carry the insert with both junctions spanned.
- **Unknown molecules:** these are mostly reads that end inside the tested span (`partial:ends`, `partial:short_read`).
  The step of each one is in `unknown_by_step`.
  - SAMN11265848 has 1,028 because its ONT reads are short at this site. The figure keeps them all.
- **No conflicts:** no molecule's records disagree, and no molecule the search saw with the core here is called A.
- **Seen only by the strict test:** 1–6 molecules per site are strict P but were not counted by the search at
  this position. Why was not traced.
- **Empty-site counts:** the v0.2 rule gives 197 / 53 / 78 / 22 molecules against v0.1's 177 / 53 / 48 / 16. v0.1
  rejected reads whose whole alignment had a net length difference over 200 bp anywhere; v0.2 checks only the site.
- **Wording:** at each of the 4 sites, some molecules of the sequenced DNA carry an insert of 1.2–1.4 kb with both
  junctions spanned by the read check's rules, and others span the empty site. These are candidate structures in
  part of the sequenced DNA. They are not insertion events, mobility, or new sites in a population sense.
- **Search class:** the search no longer classes them as `not_in_assembly:core_absent`, because their inner flank
  ends are not aligned. They stand as strict-test results at 4 sites chosen by v0.1.
- **Other such sites not looked for:** at the main setting, 8 unplaceable contexts have both flanks
  `inner_end_unaligned`.
  - With neither flank placed, these contexts are grouped by core only, not by location.
  - Two of them are the cores of the 3 ONT candidates (58 and 42 molecules).
  - The other 6 are 6 cores with 3–4 molecules each.
  - Whether other sites of this kind exist was not searched for systematically.

### Recovery misses, v0.2 (`read_site_search_pilot/recovery_misses_v02.tsv`)

The 37 round-5-supported positions not recovered:

| reason | positions |
|---|---:|
| the core's own observations there have one flank ambiguous (repeat context) | 15 |
| no two-sided observation of the core there at all | 18 |
| found with fewer than 3 molecules | 2 |
| one flank unplaced | 1 |
| 2 two-sided observations, merged into a junction and a shifted position | 1 |

- **Cores with no two-sided observation:** these are short (183–507 bp) one-gene cores lying inside a larger
  multi-gene occurrence. The larger occurrence is recovered as its own core.
- **v0.1's explanation was wrong:** none of the misses is fixed by matching, so the v0.1 explanation (containment
  vs overlap) was wrong.

### Still not done

- Placing unplaced or ambiguous flanks on other assemblies of the pool.
- A systematic search for sites whose flanks stop short of the core consistently (an insert larger than the core).
  That would need its own rule, fixed in advance.
- Any assembly outside the frozen 11; no downloads.
