# Graph entry on frozen Comamonadaceae subsets (2026-10-27)

Gap-independent: PanGraph repeated blocks from verified FNA only. No gap list, empty site, FASTQ or new
download was used. Evidence is `assembly_only`. Graph blocks are not confirmed complete insertions,
groups are not insertion families, pairs are not independent events, and organism labels are not species.
`production_ready=false`, `annotation_complete=false`, `site_groups_computed=false`.

## Subset selection (step 1)

`select_subsets.py` over the stage-1 Comamonadaceae pool:
- inputs: S2 `catalog.tsv` (sha256 `6ec4ba8b…`) and skani 0.4.0 triangle edges (`c410a2ab…`, measured ANI/AF);
- eligibility: Complete Genome, not atypical, `verify.json` all checks true, GFF3 + FAA present (412 eligible);
- seed species clusters: Delftia `SAMD00046762`, Comamonas `SAMD00018645`, Paracidovorax `SAMN02598334`;
- anchor: the member with the most near (ANI ≥ 99) partners from other BioProjects, ties broken by accession;
- members: per ANI-to-anchor bin, near [99,100], intermediate [95,99), medium [80,95), two each, same genus,
  each from a new BioProject and clone group, in accession order.

The anchor rule was changed once, before any graph was built. Taking the first accession as anchor gave no
near partner for Delftia or Comamonas. Bins that can't be filled stay short; they are not back-filled:

| subset | genomes | near | intermediate | medium | shortfall |
|---|---:|---:|---:|---:|---|
| comamonas | 5 | 2 | 0 | 2 | intermediate 2 (no eligible genome at 95–99 % ANI to the anchor) |
| delftia | 7 | 2 | 2 | 2 | none |
| paracidovorax | 4 | 2 | 0 | 1 | intermediate 2, medium 1 (only one eligible genome in range) |

Manifests are `<subset>.samples.tsv`: FNA/GFF3/FAA paths and sha256, BioProject, ANI, AF. The rule and
input hashes are in `selection.json`.

## Run (step 2)

`core_context_subsets.sbatch`, job 27327980:
- lr6, 16 CPU;
- PanGraph 1.4.0 Linux binary, sha256 `462949f3…`, the same binary as the H16/H850 smoke;
- `build --verify` linear (circular paths deferred, none excluded);
- `call` defaults: 50 ≤ block < 2500 bp, 500 bp flanks, core ≥ 0.9, similar flanks ≥ 0.9, low < 0.7.

The job ran on commit `fab5ca8` plus the uncommitted provenance change (`binary_sha256`, `build_seconds`) that
is committed together with these results. Graph JSON, joint FNA and contexts.fna stay on scratch
(`$PROJ/core_context/subsets_20261027/<subset>/`).

| subset | bp | build wall / RSS | call wall / RSS | graph sha256 |
|---|---:|---|---|---|
| comamonas | 28.3 M | 1:35 / 2.4 GB | 2:42 / 0.2 GB | `cd611a99…` |
| delftia | 45.4 M | 2:58 / 3.3 GB | 13:30 / 0.4 GB | `26047570…` |
| paracidovorax | 20.7 M | 0:52 / 1.1 GB | 3:50 / 0.2 GB | `7b2915df…` |

Pairwise Biopython alignment dominates the call step. Scale-up needs pair sampling or a faster aligner first.

## Results (step 3)

Groups are repeated graph blocks with at least one instance of 50–2500 bp. All instances are kept, including
second copies in the same assembly. No pair review was truncated.

| | comamonas | delftia | paracidovorax |
|---|---:|---:|---:|
| graph blocks / repeated (any length) | 2,380 / 1,470 | 5,106 / 3,590 | 2,595 / 1,719 |
| repeated blocks with every instance ≥ 2500 bp (not called) | 554 | 931 | 629 |
| repeated block groups in range | 914 | 2,656 | 1,087 |
| groups multi-copy in some assembly | 159 | 135 | 86 |
| instances / second copies in same assembly | 3,421 / 729 | 11,319 / 574 | 3,990 / 530 |
| context pairs | 10,629 | 25,595 | 12,625 |
| similar flanks (same-context candidate) | 5,079 | 9,890 | 2,838 |
| different adjacency + low flank identity (different-context candidate) | 601 | 1,474 | 2,479 |
| context unresolved | 3,850 | 10,027 | 6,312 |
| core identity < 0.9 | 1,096 | 4,150 | 928 |
| edge or missing context | 3 | 54 | 68 |
| groups with ≥ 1 different-context candidate pair (multi-copy among them) | 132 (64) | 257 (53) | 128 (32) |

Different-context candidate pairs by copy structure:

| subset | multi-copy group, cross-sample | multi-copy group, same sample | single-copy-per-genome group, cross-sample |
|---|---:|---:|---:|
| comamonas | 232 | 242 | 127 |
| delftia | 336 | 336 | 802 |
| paracidovorax | 1,465 | 771 | 243 |

Full tables per subset:
- `review_summary.json`, `families.tsv`, `occurrences.tsv`;
- `context_pairs.tsv.gz`, with all three relation classes kept;
- `boundary_flags.tsv`, `case_pairs.tsv`, `graph_provenance.json`, `contigs.tsv`.

## Review figures

Each subset has a `cases.svg`. It shows one similar, one different-candidate and one unresolved pair:
- coloured 500 bp U/D, core in orange;
- ±4 kb GFF3 CDS (blue), RNA (red) and pseudo/partial (grey) drawn for review only;
- block direction unified.

The pairs are chosen by a fixed rule, not by eye: the group with the most samples, a cross-sample pair if one
exists, then the smallest ids. They are therefore typical of the biggest groups, not the best candidates:
- **comamonas:** all three cases are rRNA-operon blocks (~1,160 and 1,560 bp).
- **delftia:** the different-context case is EF-Tu (`tuf`), one complete CDS in each copy. These are the known
  tufA/tufB paralogs; this is a paralog pair, not an insertion.
- **paracidovorax:** the different-context case is an IS5-family transposase block (752 bp). One copy contains the
  complete CDS; in the other copy a core end cuts a feature.

So low flank identity plus different adjacency also catches rRNA operons and housekeeping paralogs.
Separating these needs the protein/annotation layer and boundary review.

## Boundary problems (step 4)

These are GFF3 overlap previews, not the S4 coding check:

| instances | comamonas | delftia | paracidovorax |
|---|---:|---:|---:|
| next to another repeated block | 3,149 (92 %) | 10,920 (96 %) | 3,721 (93 %) |
| a core end cuts a CDS/RNA (> 20 bp from its ends) | 2,229 | 8,043 | 2,983 |
| 1–2 complete CDS inside, no cut | 506 | 1,888 | 554 |
| no feature overlap | 306 | 811 | 226 |
| edge-limited (< 500 bp flank) | 2 | 10 | 5 |

- **Fragmented blocks.** PanGraph splits both the shared backbone and multi-copy elements into many adjacent
  blocks. Most "repeated" groups are pieces of single-copy orthologous sequence, and a mobile element can span
  several blocks. A block boundary is therefore not an element boundary, and the different-adjacency criterion
  partly reflects how graph fragmentation varies between copies.
- **Length range.** Repeated blocks ≥ 2500 bp are not called at all (554 / 931 / 629). Larger elements
  appear only as fragments or not at all.
- **Next step for boundaries.** Merge runs of adjacent blocks with co-occurring copy number into candidate cores
  before comparing contexts, then use CDS boundaries (S4 checks). The protein entry is independent of graph
  boundaries and is the next step.
