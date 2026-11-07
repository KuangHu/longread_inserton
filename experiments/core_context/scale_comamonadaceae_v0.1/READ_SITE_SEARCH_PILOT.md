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

### Still not done (at v0.2)

- Placing unplaced or ambiguous flanks on other assemblies of the pool.
- A systematic search for sites whose flanks stop short of the core consistently (an insert larger than the core).
  That would need its own rule, fixed in advance. **Done in v0.3 (below).**
- Any assembly outside the frozen 11; no downloads.

## v0.3 (2026-11-02): the core hit is a seed, not the insertion boundary

Requested after the 4 v0.1 sites were found only by hand. The fix is a general rule, not an allowance for those
sites; thresholds are unchanged. Rules are frozen at git tag `read_site_search_v0.3_frozen` (commit 650b5c8), before
the validation batch was chosen.

### The rule (`boundary`, `flank_of`, `insert_fields` in `../read_site_search.py`)

1. **Search window:** per core hit and side, the read beyond the hit, up to 2,500 − hit length + 2 kb, is aligned to
   the sample's own assembly. The 2,500 bp bound is the project's target insertion size (`match_both_lt2500`)
   applied to the whole putative insertion, not a per-case allowance.
2. **Host anchor:** the alignment nearest the core whose aligned read span is at least 80% of L. Its inner end is
   the putative junction.
   - The read bases between the core hit and the junction are **extra sequence**, recorded per side. Up to
     20 bases counts as 0 (alignment end noise, `INNER_CLIP`).
   - Shorter alignments in between stay extra. An example is an element end that also occurs elsewhere.
   - The search does not skip an inner-most anchor that sits in a repeat in favour of a farther unique alignment,
     because that would move the junction. Such a position becomes unplaceable instead.
3. **Placement:** exactly L bp of read ending at the junction is the flank, placed by the v0.2 rule. The host anchor
   must therefore itself cover the junction: an unaligned inner end is still `inner_end_unaligned`, never ignored.
4. **Per merged position:**
   - **Lengths reported:** core hit length, extra sequence per side (median and range, assembly orientation), and
     the putative insertion length (`insert_len`, median of extra + core hit + extra).
   - **`boundary`:**
     - `unknown` unless, on each side, at least 80% of the molecules have an extra length within 50 bp of that
       side's median;
     - otherwise `exact` when both medians are 0, else `consistent`.
   - **`insert_size`:** `le_2500` or `gt_2500` of the putative insertion (since v0.3.2: `lt_2500`, strictly
     < 2,500 bp, or `ge_2500`). Whether it holds 1–2 genes is not checked here.
   - **The insertion replaces the core:** its length replaces the core length in the geometry, class and empty-site
     rules.
5. **Empty site:** still auxiliary evidence, never a discovery condition.

**Changed on the development samples before freezing** (each with a test in `../test_read_site_search.py`):
- **Boundary rule:**
  - First version: all molecules within 50 bp.
  - Problem: 106 of the known positions came out `unknown` because of one stray molecule each, though their
    median extra was 0.
  - Final: the 80% majority above.
- **Junction from two strands:** a plus/minus junction read from its two strands was two positions. It is now
  merged (`swapped`, `canonical`).

**Counterexample tests:**
- extra sequence beyond the core;
- partial core hit;
- a short alignment between anchor and core;
- repeated flank (ambiguous, not skipped);
- neighbouring different insertion, both absent from and present in the assembly;
- one stray molecule;
- plus/minus junction from both strands.

### Regression on the 11 development samples (flanks 1 kb, ≥ 3 molecules), v0.2 → v0.3

| class | v0.2 | v0.3 |
|---|---:|---:|
| `known_catalog_position` | 530 (all `exact` at v0.3) | 530 |
| `known_catalog_position_other_core` | 600 | 601 (600 `exact`) |
| `catalog_mismatch:coverage` | 64 | 65 |
| `catalog_omission:no_annotated_cds` | 6 | 6 |
| `homolog_not_in_catalog` (group / other) | 138 / 73 | 135 / 75 |
| `not_in_assembly:core_absent` | 0 | **5** (all `consistent`, `le_2500`) |
| `not_in_assembly:junction` | 0 | 1 |
| `contig_end_join` | 3 | 8 |
| `unplaceable_context:*` | 390 | 109 |
| `shifted_known_position` | 92 | 98 |

- **Recovery:** round-5-supported 329 / 366, others 201 / 227 (both unchanged).
- **Negative control:** 47 genes with exactly one position, 8 with none, none with two or more (unchanged). The
  8 control-gene contexts that were unplaceable at v0.2 are gone.
- **Unplaceable contexts:** 390 → 109. The 243 v0.2 contexts with an unaligned inner flank end no longer exist as a
  class: their junctions are now found beyond the extra sequence. Of the 109 left, 99 have a flank in a repeat
  (`ambiguous`) and 10 have flanks without an anchor or unplaced.
- **The 4 development cases are found automatically,** each `consistent`:

  | site | molecules | extra L / R (median, range) | insert / core hit |
  |---|---:|---|---|
  | SAMN30525180 3,553,051–059 | 41 | 69 (58–156) / 36 (30–54) | 1,326 / 1,220 |
  | SAMN30525180 1,653,143–151 | 9 | 69 (65–75) / 35 (33–39) | 1,323 / 1,217 |
  | SAMN11265848 3,064,452–457 | 23 | 21 (0–60) / 76 (68–94) | 1,176 / 1,078 |
  | SAMN15738716 3,158,294–304 (HiFi) | 10 | 57 (56–58) / 0 | 1,395 / 1,339 |

  These 4 are development regression cases. They do not validate the rule.
- **A fifth site (also development data):** SAMN30525180 NZ_JAOCIM010000001.1:1,564,844–853, same core
  PF_a5c0db98e215 with the same extra lengths (69 / 34 bp), 5 molecules.
  - Strict P/A (`candidates/`): 7 P matched, 57 A matched, no conflict, 51 unknown.
  - The site lies inside a catalog occurrence of an unrelated core (`catalog_unmatched`).
  - The strict P/A counts of the 4 earlier sites are unchanged.
- **New junction:** SAMN15738716 PF_fd55e20062f5, 62 molecules. It joins NZ_JARXUY010000001.1:2,231,535 (+) to
  4,097,738 (−), 1.87 Mb apart in opposite orientation, through an insertion of 1,333 bp (core hit 900 + 433 bp
  extra).
  - This is a candidate structure the assembly does not show. It could be an inversion in part of the culture or an
    assembly join error.
  - It is not resolved here.
- **Contig end joins (5 new):**
  - Mostly SAMEA115676957 (contig-level) with 0–206 bp of extra sequence.
  - SAMN52101483 NZ_CM129821.1, 89 molecules, 414 bp extra.
  - SAMN15738716 contig 2, 48 molecules, 293 bp extra, the same core as its candidate.
  - These were unplaceable at v0.2 because their inner flank ends did not align.

## Validation batch 1 (2026-11-02): frozen rules, samples not used in development

- **Samples (`read_site_search_validation1/`, frozen before the run):** the pilot freeze rule, unchanged, applied
  after excluding every sample used in development (`freeze.json`: the 11 pilot samples, SAMN31258398 and
  SAMN09743972).
  - **Result:** 7 samples from 6 genera, 5 ONT and 2 HiFi: SAMN51281039 (Variovorax), SAMN30525383 (Delftia),
    SAMN38562701 (Comamonas, HiFi), SAMEA115676833 (Acidovorax), SAMN56368990 (Hydrogenophaga), SAMN43526786
    (Polaromonas), SAMN37751820 (Variovorax, HiFi).
  - **Inputs:** reads already on disk, no downloads. The panel was rebuilt for these samples' catalog and control
    genes (same 1,083 cores).
- **Code:** the search ran at the frozen rules (commit 0293d34 = the tag plus the frozen sample list).
  `candidate_evidence.py` was changed afterwards only to follow a read-check fix, with no effect on ONT/HiFi. No
  validation candidate reached it.

### Results (flanks 1 kb, ≥ 3 molecules)

| class | positions |
|---|---:|
| `known_catalog_position` | 314 (all `exact`) |
| `known_catalog_position_other_core` | 385 (all `exact`) |
| `catalog_mismatch:coverage` | 35 |
| `catalog_omission:no_annotated_cds` | 1 |
| `homolog_not_in_catalog` (group / other) | 112 / 62 |
| `not_in_assembly:*` | 0 |
| `contig_end_join` | 0 |
| `unplaceable_context:*` | 11 |
| `shifted_known_position` | 75 |

- **Recovery:**
  - round-5-supported: 181 / 190 (95%);
  - others: 133 / 135;
  - control genes as catalog entries: 17 / 17.
- **Negative control:** 28 genes searched. 25 have exactly one position, 3 have none, and none has two or more at
  every N ≥ 2. At N = 1 with 2 kb flanks, one gene has 2.
- **Candidate structures:** none reach 3 molecules.
  - Below threshold: 3 `core_absent` (1 molecule each) and 27 junctions (1 molecule each).
  - None of these comes close to the N = 3 bar (none has 2 molecules).
- **Sensitivity grid:** in `summary.json`.
  - Junctions appear only at N = 1 (26–35).
  - Recovery is 176–183 / 190 across flank lengths and N.

### Reading

- **Rule held:** on samples not used to build it, the frozen rule keeps known-site recovery (95%) and the negative
  control clean.
- **Boundaries:** known positions come out `exact` (v0.3.1: `no_extension_detected`). **Corrected in v0.3.1:**
  this is not a check of boundary accuracy. It only says that the reads join U and D next to the core hit, with ≤
  20 bp counted as 0. Where the assembly holds the whole element, element ends beyond the core align as flank.
- **No extended-insertion structure:** the batch has no structure the assembly lacks with an insertion longer
  than its core.
- **What this does not show:** this batch does not measure the rule's sensitivity for such structures. There is no
  positive case in it to find. One small batch without findings is not evidence that such structures are rare.


## v0.3.1 (2026-11-02): boundary counted by molecule; full-insertion check of the 5 candidates

### Boundary statistics (`insert_fields`)

- **What was wrong (found in review, reproduced with the function itself):**
  - v0.3 took medians and the 80% agreement over observation records, not molecules.
  - It also tested each side separately.
  - Counterexamples, now tests:
    - one molecule with 8 identical records plus 2 molecules with different ones was called `exact`, with only
      1 / 3 molecules supporting both sides;
    - 80% agreement on each side but 60% on both at once was also called `exact`.
- **Fix:**
  - Each molecule's records are reduced to one value per side (their median) before the medians over molecules.
  - A molecule agrees on a side if its extra length is within 50 bp of that median.
  - The boundary needs at least 80% of the molecules to agree on both sides at once.
  - Reported per position: `boundary_molecules`, `left_agree_molecules`, `right_agree_molecules`,
    `both_agree_molecules`.
- **Labels:**
  - `exact` → `no_extension_detected`;
  - `consistent` → `extension_consistent`;
  - `unknown` is unchanged.

  `no_extension_detected` means only that the reads join U and D directly next to the core hit; ≤ 20 bp of extra
  sequence counts as 0. It does not establish the element's own boundaries: where the assembly holds the whole
  element, element ends beyond the core align as flank.
- **Effect on real data:** none. Rerun on both batches; jobs are in each `jobs.txt`.
  - Every position keeps its boundary class: development 530 / 600 known positions `no_extension_detected`;
    validation 314 / 385.
  - The 5 candidates:

    | site | molecules | agree on both sides |
    |---|---:|---:|
    | SAMN30525180 3,553,051 | 41 | 40 |
    | SAMN30525180 1,653,143 | 9 | 9 |
    | SAMN30525180 1,564,844 | 5 | 5 |
    | SAMN11265848 | 23 | 23 |
    | SAMN15738716 | 10 | 10 |

    All 5 are `extension_consistent`.

### Full-insertion check (`../candidate_structure.py`; `read_site_search_pilot/candidates/structure_*`)

Three levels, reported apart:
- **core hit:** coverage, identity and both ends of the panel core sequence in the insertion;
- **structure:** molecules joining U – insertion – D (strict P/A) and the boundary counts by molecule;
- **complete target:** the insertion between the junctions is < 2,500 bp and holds 1–2 complete genes.

Method:
- **Insertion sequence:** a consensus from every record of the molecules the strict check called P, never one read.
  - Draft: the P allele window, the insertion ± 300 bp of assembly.
  - Polishing: minimap2 + `samtools consensus` (simple mode, call fraction 0.6, min depth 3), 2 rounds.
  - Junctions: blastn of 300 bp of U and D against the consensus. A target-site duplication stays in U and D.
- **Genes:**
  - ORFs on both strands (table 11, ≥ 50 aa), blastp against the sample's annotated proteins and every pool protein
    of the core's family.
  - A complete gene is a complete ORF whose best hit covers ≥ 90% of the subject and of the ORF.
  - Ambiguous consensus bases (N) are counted. If any lies outside a complete gene, the gene count is not
    interpreted (`unresolved`). A frameshift in a low-depth consensus is never read as a broken gene.
- **Assembled copies:** the insertion is aligned to the sample's own assembly (asm20, all copies). Copies covering
  ≥ 90% of it are listed with the annotated CDS inside and their catalog status.

| site | P molecules (records) | insertion bp | depth min / median | core coverage / identity / both ends | complete genes | gene (best hit) | genes status | assembled copies (in catalog) | target (< 2,500 bp, 1–2 genes) |
|---|---|---:|---|---|---:|---|---|---|---|
| SAMN30525180 3,553,051 | 46 (46) | 1,332 | 37 / 46 | 100% / 99.84% / yes | 1 | IS256 family transposase, 408 aa, 99.8% | 2 ambiguous bases, inside the gene | 16 (15) | met |
| SAMN30525180 1,653,143 | 10 (10) | 1,331 | 8 / 12 | 100% / 99.84% / yes | 1 | IS256 family transposase, 408 aa, 99.8% | 2 ambiguous bases, inside the gene | 16 (15) | met |
| SAMN30525180 1,564,844 | 7 (7) | 1,328 | 5 / 7 | 100% / 98.70% / yes | 0 | (a 325-aa ORF covering 80% of the IS256 protein) | **unresolved**: 19 ambiguous bases, 1 inserted base against the assembled copy | 16 (15) | unresolved |
| SAMN11265848 | 30 (30) | 1,177 | 25 / 32 | 100% / 100% / yes | 1 | IS630 family transposase, 359 aa, 100% | resolved | 23 (21) | met |
| SAMN15738716 (HiFi) | 12 (12) | 1,402 | 11 / 12 | 100% / 100% / yes | 1 | IS4 family transposase, 445 aa, 100% | resolved | 12 (12) | met |

- **Core level:** every insertion holds the whole panel core, both ends covered.
- **Gene count:** the gene count comes from the insertion, not from the seed core. Each resolved insertion holds
  exactly one complete gene, the core's own transposase family.
- **The 1,564,844 site:**
  - At 5–7 molecules the consensus is not good enough to call genes. Its one inserted base against the assembled
    copy and its 19 ambiguous bases explain the shortened ORF; the gene is not called broken.
  - Its sequence otherwise matches the other SAMN30525180 insertions and 16 assembled copies at 99.9%.
- **Ambiguous positions 252–253** at the other two SAMN30525180 sites: the molecules there carry mixed bases at
  37–46× depth. This could be a two-base difference among the copies or a systematic read error. It is not
  resolved here.
- **Assembled copies:** every insertion matches 12–23 assembled copies elsewhere in its own genome at 99.9–100%.
  Nearly all of them are catalog occurrences of the same core.
- **Reading:**
  - At these 5 sites, part of the sequenced DNA carries a copy of a multi-copy insertion-sequence-family element at a
    site where the assembly shows the empty site.
  - Four meet the complete target (< 2,500 bp, one complete gene); one is unresolved.
  - This is a candidate structure in part of the DNA. It is not an independent insertion event, and mobility is not
    shown.
- **Not an insertion:** the 62-molecule junction of SAMN15738716 (PF_fd55e20062f5, 2,231,535 + ↔ 4,097,738 −, an
  insertion of 1,333 bp between) is listed apart as a structural anomaly (`candidates/structural_anomalies.tsv`).
  It is not treated as an insertion.

### v0.3.2: target size strictly < 2,500 bp

- **The change:** the size class of a putative insertion is now `lt_2500` (< 2,500 bp) or `ge_2500`; before it was
  `le_2500` / `gt_2500`, so exactly 2,500 bp counted as in target. The complete-target column of the structure
  check uses the same rule (`size_class`).
- **Boundary test:** 2,499 bp → `lt_2500`; 2,500 and 2,501 bp → `ge_2500`, for the insertion length (extra + core
  hit + extra), not the core hit.
- **Effect on earlier results:** no merged position in the development or validation-1 runs has an insertion of
  exactly 2,500 bp, so only the label changes. No other threshold changed.

### Next (decided 2026-11-02: proceed; see the scale-up section)

- **Second validation batch:** after these rules are frozen, as a new tag.
- **Computational positive test,** fixed in advance and reported apart from real samples:
  - in validation samples, cut a seeded random set of known catalog occurrences out of the assembly, with an extra
    0–300 bp drawn per side, and search the real reads against the edited assembly;
  - the cut sequence is the truth for the site and for the insertion length and boundary;
  - this measures sensitivity for insertions longer than the core, which a random batch without positives cannot.

## Scale-up (2026-11-02 to 11-06): all eligible Comamonadaceae samples, frozen v0.3.2

### Freeze

- **Code:** tag `read_site_search_v0.3.2_frozen` (commit 3b72e64). Relative to v0.3.1 the only rule change is that
  the target size is strictly < 2,500 bp. It is its own commit, with boundary tests at 2,499, 2,500 and 2,501 bp.
  - No threshold was changed after any scale result.
  - Later commits only add bookkeeping: site layer, acceptance, cleanup, positive test, retry-aware acceptance.
- **Per-batch record:** each batch's `freeze.json` stores:
  - the eligibility rule text and the parameters;
  - the code commit;
  - the ready-manifest path and its sha256 (`ba60fe83b9db…` for every batch);
  - the sample list with run sha256s.
- **Main setting:** flanks 1 kb, at least 3 independent molecules. The core hit is a seed; no empty site is required;
  a single read is never a site.
- **Batches:**
  - validation 2 (6 samples) and scale 1–3 (40 each) are every eligible sample not used in development, with
    genera interleaved round-robin;
  - development (11) and validation 1 (7) were rerun on v0.3.2: labels changed, no position, class, count or length
    did.

### Downloads

`download_batches/site_search_download_check.py` → `site_search_download_check.tsv`: **no run can be added to this
search.** Every usable ONT/HiFi run of a manifest sample is already verified. The rest are:

| reason | runs | BioSamples | size |
|---|---|---|---|
| no assembly in the reviewed manifest (nothing to place flanks on) | 113 | 101 | 97.4 GB |
| assembly excluded as NCBI atypical | 3 | 3 | 21.1 GB |
| run not usable (PCR / amplicon / WGA / single cell / metagenome / mean read < 500 bp / no ENA md5) | 49 | 47 | 22.6 GB |

So nothing was downloaded in this round, and compute read only the verified ready manifest. Of the 231 samples
considered, 87 were not eligible and are not in the search (`read_site_search_scale3/candidates_considered.tsv`):
- 49 with CLR or reads not countable as molecules, alone or with another reason;
- 33 below 30× depth;
- 7 with other read types;
- 2 aligned-BAM submissions.

The reasons overlap; see the file for exact combinations. CLR and possibly trimmed data were not mixed in.

### Results by batch

Generated tables: `read_site_search_scale_report.md` (`../scale_report.py`). Real batches, development and the
computational test are never summed together. All batches are complete by sample ID: 0 failed, 0 awaiting input.

**Real samples not used in development (validation 1–2, scale 1–3; 133 samples):**

| | validation1 | validation2 | scale1 | scale2 | scale3 | total |
|---|---|---|---|---|---|---|
| samples complete | 7/7 | 6/6 | 40/40 | 40/40 | 40/40 | 133/133 |
| known recovery, round-5 supported | 181/190 | 151/166 | 761/808 | 1336/1405 | 757/794 | 3186/3363 (94.7%) |
| known recovery, other catalog | 133/135 | 183/197 | 753/806 | 1225/1311 | 973/990 | 3267/3439 (95.0%) |
| known sites | 617 | 645 | 3397 | 5653 | 4054 | 14366 |
| in assembly, not in catalog | 149 | 94 | 809 | 1026 | 1083 | 3161 |
| **not-in-assembly candidate sites** | 0 | 0 | 8 | 1 | 1 | **10** |
| candidate positions / core keys / samples | 0 | 0 | 9 / 8 / 4 | 1 / 1 / 1 | 2 / 2 / 1 | 12 / 11 / 6 |
| structural anomaly (distant / inverted junction) | 0 | 1 | 66 | 9 | 7 | 83 |
| unresolved | 28 | 28 | 184 | 206 | 154 | 600 |
| boundary unknown (counted sites) | 1 | 3 | 24 | 2 | 18 | 48 |
| ambiguous placement (repeat) / truncated alignment | 11 / 0 | 130 / 0 | 358 / 6 | 1068 / 0 | 143 / 1 | 1710 / 7 |
| max RSS (GB) | 1.43 | 11.01 | 11.38 | 11.29 | 15.55 | |
| cleanup freed (GB) | 0.12 | 0.29 | 1.58 | 1.67 | 2.36 | 6.02 |

Notes on the table:
- Counts are physical sites, with cores at one place merged; candidate positions and cores are counted separately.
- Core keys: 11 distinct in total.
- "Unresolved" means boundary unknown, repeated placement or insufficient coverage.
- Ambiguous placement and truncated alignment are unplaceable contexts with ≥ 3 molecules.

**Negative controls (single-copy housekeeping genes), by physical site:**
- Validation 1–2 and scale 2: every searched gene at one site.
- Scale 1: 265 at one site, 1 at none, 7 at two.
  - 5 of the 7 are one locus whose clusters differ by more than 50 bp on one junction.
  - The other 2 are separate loci: SAMN51281032 RNAP β/β′ also sits on plasmid NZ_JBRAYR010000004.1, a real second
    assembled copy.
- Scale 3: 264 at one site, 1 at none, 3 at two, 1 at three. All 4 multi-site genes are
  `one_locus_junction_differs_gt_TOL`.

**The 10 real candidate sites (12 positions).** These are "candidate structures in part of the sequenced DNA", not
insertion events. Strict per-molecule P counts are 0 throughout; the P allele is only matched (P_matched), and the
empty allele is auxiliary.

| batch | site | core(s) | P_matched / A_matched | insertion | full-insertion check |
|---|---|---|---|---|---|
| scale1 | SAMEA6863324 NZ_LR813086.1:398,381 | PF_1d912f750caa | 6 / 373 | 1,029 bp | unresolved, ambiguous consensus |
| scale1 | SAMEA6863324 NZ_LR813086.1:3,619,182 | PF_31cb679e911e (+PF_0b44bc63f72c) | 3 / 347 | 835 bp | unresolved, ambiguous consensus |
| scale1 | SAMEA6863356 NZ_CAHPRW010000010.1:10,237 | PF_09076b68b6fc | 3 / 569 | 1,023 bp | unresolved, ambiguous consensus |
| scale1 | SAMEA6863356 NZ_CAHPRW010000010.1:11,355 | PF_09076b68b6fc | 3 / 572 | 1,240 bp | unresolved, ambiguous consensus; boundary unknown |
| scale1 | SAMN28084394 NZ_CP097269.1:821,263 | PF_12774aadd19f | 6 / 45 | 1,154 bp | unresolved, ambiguous consensus |
| scale1 | SAMN28084394 NZ_CP097269.1:820,188 | PF_6eeae208e0a0 | 3 / 52 | 1,261 bp | unresolved, ambiguous consensus |
| scale1 | SAMN42749333 NZ_CP162616.1:5,467,735–5,467,993 | PF_5194ba28f00a | 3 / 158 | — | junction not found in consensus |
| scale1 | SAMN42749333 NZ_CP162616.1:456,191–458,920 | PF_7dd788787bef | 4 / 85 | 2,473 bp (assembly: 2,729 bp of other sequence) | unresolved, ambiguous consensus (gene range 1–4); reviewed separately, see below |
| scale2 | SAMEA120309485 NZ_CDRKJB010000001.1:3,210,811 | PF_d6e933620b0d | 6 / 114 | 1,302 bp | unresolved, ORFs not explained by a complete gene (range 0–2) |
| scale3 | SAMN13494298 NZ_CP046622.1:3,263,988–3,264,147 | PF_052544d04561; PF_2f595d1c0409−,PF_82422e478eda+ | 21 / 577 | 2,150 bp | unresolved, ambiguous consensus |

About these candidates:
- None reaches "complete target met". Most are rare structures in 3–6 molecules against 45–577 empty-allele
  molecules, so the consensus is thin.
- The scale-2 insertion has 3 assembled copies elsewhere in the same genome. Its core is 100% identical, but no ORF
  is a complete gene.
- The scale-3 site:
  - has 21 molecules;
  - has 26 ambiguous bases in a 2,150 bp consensus;
  - its two seed cores cover only 55% and 90% of the core, at 86–88% identity, with neither end covered;
  - its gene count is unresolved (0–3).

### Computational positive test (`read_site_search_positive1`; reported apart from real samples)

**Design:**
- Seed 20261102. 100 known catalog occurrences cut from the assemblies of the 13 validation samples (`positive_test.py
  prepare`).
- Strata: core length bin × single / two-gene × repeated / single-copy × extra sequence (none, or 1–300 bp per side).
- Edited assemblies, catalog, coordinates and genome alignments are rebuilt. The truth (`read_site_search_positive1_truth/`,
  sha256 in `prepare.json`) is used only by `score`.

**Result:** 99/100 cuts detected as not_in_assembly with ≥ 3 molecules. Every length bin, structure, repeat
class and extra mode is at 100%, except one single-copy < 1 kb single-gene cut.
- **The miss** (SAMN37751820 PF_4aecb2fa43e0): it was found with 126 molecules, but at a 51 bp junction offset
  (2,413,797 / 2,413,848). That is outside the 50 bp tolerance, so it stays scored as not detected.
- **Insertion length:** 72 cuts within 10 bp, 27 within 11–50 bp. The size class agrees in 99/99.
- **Extra sequence per side:** 92 within 20 bp, 2 within 21–50 bp, 5 above 50 bp. In all 5, the total insertion
  length is still within 20 bp. The left/right split is measured from whichever core seeded the read; these are
  two-gene or repeated windows.
- **Boundary:**
  - 47 extension_consistent;
  - 50 no_extension_detected (no evidence of extension; this is not boundary validation);
  - 2 unknown.
  - Against the truth: cuts with extra sequence 45 extension_consistent / 1 no_extension_detected / 2 unknown (of
    48); cuts without extra sequence 49 no_extension_detected / 2 extension_consistent (of 51).
- **Consensus:** 49 met, 25 unresolved with ambiguous consensus, 25 unresolved with unexplained ORFs.
  - Consensus length is within 10 bp of the truth for 98.
  - Unresolved here is the rule doing its job (an unexplained ORF keeps the gene count open), not a miss.
- **Not-in-assembly calls away from cuts:** 24.
  - 9 are also in the unedited real run;
  - 15 have one end within 2.5 kb of a cut and are side effects of the edit;
  - 0 unexplained (`positive_test_away_from_cuts.tsv`).

### Failures, retries and integrity

- **Scale 3:** 5 of 40 search tasks failed in the first run (tasks 19, 25, 27, 31, 38).
  - All failures ran on n0334.lr6.
  - The error was Errno 108, "cannot send after transport endpoint shutdown", while the ready-manifest check stat'ed
    read files on scratch. The files are intact and readable from login.
- **Retry 1** (n0334 excluded): 19, 25 and 31 completed; 27 and 38 failed on n0348 / n0357.
- **Retry 2:** 27 and 38 failed on n0258 / n0194.
- **Retry 3** (30 min later, 5 nodes excluded): both completed. Combine completed. Evidence failed on n0184 with the
  same error, and its rerun completed, as did structure and sites.
- The 35 search tasks that succeeded first time were not rerun. Method and code were unchanged.
- Every attempt is in `read_site_search_scale3/jobs.txt`. `acceptance.json` lists each failed or cancelled attempt
  under `retried`, with the job that replaced it.
- Acceptance counts an attempt as retried only if a later job of the same stage and array index completed; otherwise
  it stays a failure (test in `test_scale_sites.py`).
- Final state: 40/40 complete, `failed_jobs` empty.
- The other batches ran without failures.

### Cleanup

After acceptance, each batch's work directory is archived (`archive.tar`, sha256 in `archive.json`). Then the
rebuildable files listed in `cleanup_manifest.tsv` are deleted. They are listed by explicit pattern under the
absolute work path:
- boundary / placement / panel alignments;
- consensus rounds and drafts;
- reads extracted from the verified FASTQ;
- blast temporaries.

Kept:
- per-molecule evidence (`positions`, `observations`, `search.json`);
- the panel;
- candidate allele windows and alignments;
- consensus `final.fa` / `final.bam`;
- the edited assemblies and their genome alignments of the positive test.

Read files, manifests and the repository are never touched. Space reclaimed:
- real batches: 6.02 GB;
- development rerun: 0.42 GB;
- positive test: 0.73 GB;
- total: 7.17 GB.

### Unresolved

- All 10 real candidate sites stay unresolved at the full-insertion level, mostly because the consensus is thin.
  More molecules would need more sequencing of those strains, not a looser rule.
- **The one `other_sequence` candidate** was reviewed afterwards, with unchanged thresholds, in
  `read_site_search_scale1_review_other_sequence/` (`candidate_evidence.py --classes not_in_assembly:other_sequence`,
  commit 5deda86).
  - Molecules: 4 P_matched and 85 A_matched; strict P is 0.
  - Boundary: extension_consistent in 3/3 molecules.
  - Insertion: 2,473 bp (lt_2500) where the assembly holds 2,729 bp.
  - Core: fully covered, 99.4% identity, both ends.
  - Consensus: depth 2–4 with 9 ambiguous bases, so the gene count is unresolved (1–4).
  - The batch run itself only reviews `core_absent` by default.
- **One locus can split into two clusters:** some negative-control loci split when one junction differs by more
  than 50 bp. They are labelled `one_locus_junction_differs_gt_TOL`, not merged. The single positive-test miss is
  the same 51 bp effect.
- 1,710 contexts with ≥ 3 molecules in real batches are ambiguously placed in repeats, and 600 sites are unresolved.
  These are reported, not called absent.

### Remaining data

- **Downloadable for this search:** 0 runs.
- **Not eligible (87 samples):** CLR or not molecule-countable, < 30×, other read types, aligned BAM.
- **Runs without a usable assembly:** 116 runs (104 BioSamples).
- **Unusable runs:** 49 runs (47 BioSamples).
- These need either an assembly (not this method's scope) or a separate decision on CLR.
