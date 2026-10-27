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
- **Next step for boundaries.** Equal copy number is not enough to merge adjacent blocks. A run of blocks becomes
  one candidate core only if the corresponding instances keep the same adjacency, order and orientation across
  every copy. CDS boundaries (S4 checks) test whether the coding units inside are complete; they are not insertion
  boundaries. Whether a context also occurs without the core is auxiliary evidence only, never an admission
  criterion (no empty site is required).

## Protein entry v0.2 (all FAA/GFF3; `protein_context.py`, job 27328730, commit 18d2d69)

v0.1 (job 27328469, `c772861`) had three faults, reproduced with counterexamples in review. They are fixed in
`9a57f44`, each with a test; its counts are superseded:
1. **Complete CDS ignored coding flags.** A CDS whose translation matched the FAA still counted as complete with
   `no_start`. Now the S4 rule applies: not partial, not pseudo, not origin-crossing, translation matches, and no
   coding flag. In these RefSeq subsets the flag rule excludes 0 CDS; the origin rule excludes 2 / 5 / 1.
2. **Unknown neighbours counted as different.** A side whose neighbours had no homology label could still yield
   `different_context_candidate`. Now a side without a shared label must be fully labelled in both occurrences;
   otherwise the pair is `insufficient_context_labels`.
3. **No length target.** Cores of 3,000 bp stayed in. Each occurrence now has a span (first start to last end; a
   two-gene span includes the spacing). Everything is computed twice: scope `all` (wide protein search) and scope
   `target_lt2500` (span < 2500 bp; the 1–2 complete coding unit target).

How the entry works:
- **Start:** every GFF3 CDS and RNA, with no graph blocks, gap list or reads.
- **Coding check:** the S4 coding check on every CDS.
- **Clustering:** MMseqs2 17.b804f at the frozen S4 settings (family 0.9/0.8, homology group 0.5/0.8), with every
  member aligned to its final representative. 0 failures, so no repairs.
- **Cores:** single-gene = one complete CDS, key = protein family; two-gene = two consecutive complete CDS, key =
  families in order with relative orientation (canonical over strands), spacing per occurrence.
- **Context:** K=3 neighbours per side in core orientation, labelled by homology group. Windows wrap only on
  `Is_circular=true` contigs.
- **Cost:** < 40 s per subset.

| | comamonas | delftia | paracidovorax |
|---|---:|---:|---:|
| genes (CDS + RNA) / complete CDS (S4 rule) | 26,338 / 25,285 | 41,878 / 39,793 | 18,892 / 18,319 |
| distinct proteins / families / homology groups | 19,420 / 11,177 / 8,005 | 26,523 / 13,101 / 9,689 | 11,533 / 6,593 / 5,892 |

Single-gene cores (wide search `all` → target `<2500 bp`):

| | comamonas | delftia | paracidovorax |
|---|---:|---:|---:|
| core keys with ≥ 2 occurrences | 4,939 → 4,816 | 6,659 → 6,449 | 4,334 → 4,189 |
| occurrences | 19,048 → 18,560 | 33,352 → 32,353 | 16,060 → 15,512 |
| multi-copy in some assembly | 231 → 224 | 195 → 182 | 86 → 79 |
| pairs similar | 26,910 → 26,207 | 69,751 → 67,786 | 20,968 → 20,256 |
| pairs different-context candidate | 600 → 599 | 777 → 756 | 1,967 → 1,956 |
| pairs one side shared | 1,258 → 1,215 | 3,771 → 3,644 | 718 → 678 |
| pairs insufficient labels | 825 → 803 | 2,705 → 2,642 | 2,031 → 2,007 |
| **core keys with ≥ 1 different-context candidate** | **52 → 51** | **85 → 79** | **24 → 21** |

Two-gene cores (`all` → `<2500 bp`):

| | comamonas | delftia | paracidovorax |
|---|---:|---:|---:|
| core keys with ≥ 2 occurrences | 4,708 → 3,503 | 6,739 → 5,022 | 4,244 → 3,087 |
| pairs different-context candidate | 65 → 58 | 110 → 99 | 27 → 17 |
| pairs insufficient labels | 372 → 322 | 1,450 → 1,138 | 237 → 136 |
| **core keys with ≥ 1 different-context candidate** | **13 → 10** | **30 → 27** | **9 → 6** |

There were no edge-limited or truncated pair reviews. Compared with v0.1, the single-gene cores with a
different-context candidate in scope `all` fell from 84 / 142 / 29 to 52 / 85 / 24, mostly because of rule 2.

Product keywords for target-scope single-gene cores with a different-context candidate are a reading aid only. "Other"
includes the intron-RT and translation/ribosomal hints; see `protein_review.json`:

| | transposase/IS | integrase etc. | hypothetical | other |
|---|---:|---:|---:|---:|
| comamonas | 13 | 3 | 4 | 31 |
| delftia | 10 | 0 | 13 | 56 |
| paracidovorax | 8 | 1 | 3 | 9 |

Overlap with the graph entry (target-scope single-gene cores with a different-context candidate; does any
occurrence lie in a graph group with a graph different-context pair?):

| | IS hint: yes / no | other: yes / no |
|---|---:|---:|
| comamonas | 11 / 2 | 13 / 25 |
| delftia | 9 / 1 | 31 / 38 |
| paracidovorax | 8 / 0 | 9 / 4 |

### U–core–D case review (target scope; `case_occurrences.tsv`, `cases_protein.svg`)

The fixed rule takes the top three single-gene cores by number of context classes linked by different-context
pairs, then one occurrence per class.

| subset | core | span | occurrences / samples / classes |
|---|---|---:|---|
| comamonas | PF_147add17340c ISL3 transposase | 1,161–1,221 | 29 / 3 / 29 |
| comamonas | PF_4a2352cf6318 IS481 transposase | 966–987 | 6 / 2 / 6 |
| comamonas | PF_2d92c93a9d2c hypothetical protein | 255 | 8 / 4 / 7 |
| delftia | PF_3070024975f6 IS481 transposase | 1,188 | 21 / 5 / 13 |
| delftia | PF_4b70bb8e2b9d IS5 transposase | 809–811 | 30 / 4 / 18 |
| delftia | PF_ad056eba1256 IS5 transposase | 793 | 15 / 1 / 14 |
| paracidovorax | PF_34809b3c6b4f ISAav5 (IS481-like) | 948 | 54 / 3 / 36 |
| paracidovorax | PF_2bf098df0d79 IS3 transposase | 1,565 | 55 / 3 / 43 |
| paracidovorax | PF_f094f2564f99 IS5 transposase | 813 | 30 / 4 / 22 |

Reading of the cases:
- **The intended pattern.** Eight of nine are IS transposase CDS: the same protein core (identical span within a
  family) between unrelated chromosomal genes. Examples are flagellar genes vs nitrogen regulation vs UvrA
  neighbours (ISL3, Comamonas), and IclR / penicillin-binding protein / tryptophan synthase neighbours (IS481, Delftia).
- **Not all are cross-strain.** Many different-context pairs are copies inside one assembly. PF_ad056eba1256 has
  15 copies, all in SAMN15905291. That is "same core, different context", but not evidence of polymorphism between
  strains.
- **Mobile neighbours.** In Paracidovorax, IS3 / IS5 copies often sit next to other mobile genes (IS5, ISAav4, Tn3,
  relaxase / mobilisation, RHS), and one IS5 copy sits next to CRISPR-cas genes. Such a context may be another mobile
  element or a nested insertion, not a host site.
- **Core is only the CDS.** The core is the transposase CDS only; IS inverted repeats and flanking non-coding parts
  are outside it. A 1,565 bp IS3 CDS suggests a fused orfA/orfB annotation; this is not checked here.
- **Unresolved hypothetical.** The 255 bp hypothetical protein (Comamonas) sits next to an AlpA phage regulator,
  arsenic-resistance genes and an IS5 copy. It is not identifiable as a mobile element from annotation and stays
  unresolved.

Outputs per subset in `<subset>/protein_entry/`:
- `protein_summary.json`, `protein_cores.tsv` (both scopes), `protein_occurrences.tsv.gz` (span, target flag,
  oriented neighbour labels / gene ids, context class per scope), `protein_pairs.tsv.gz` (all relations, both scopes);
- `member_check_repairs.tsv` (empty);
- `different_context_cores.tsv`, `protein_review.json`, `case_occurrences.tsv`, `cases_protein.svg`.

Still open:
- boundary review beyond the CDS (adjacent-block merging under consistent adjacency, order and orientation);
- paralog / rearrangement / nested-element interpretation;
- separating same-assembly copies from cross-strain contexts in the counts;
- site groups (S6);
- read checks.

Evidence stays `assembly_only`, `production_ready=false`.
