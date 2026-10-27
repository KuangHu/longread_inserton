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
| delftia | PF_3070024975f6 IS481 transposase | 1,164–1,188 | 21 / 5 / 13 |
| delftia | PF_4b70bb8e2b9d IS5 transposase | 782–811 | 30 / 4 / 18 |
| delftia | PF_ad056eba1256 IS5 transposase | 793 | 15 / 1 / 14 |
| paracidovorax | PF_34809b3c6b4f ISAav5 (IS481-like) | 948 | 54 / 3 / 36 |
| paracidovorax | PF_2bf098df0d79 IS3 transposase | 1,565 | 55 / 3 / 43 |
| paracidovorax | PF_f094f2564f99 IS5 transposase | 813 | 30 / 4 / 22 |

Reading of the cases:
- **The intended pattern.** Eight of nine are IS transposase CDS: the same protein core (one protein family; spans
  differ by at most 60 bp within a family) between unrelated chromosomal genes. Examples are flagellar genes vs nitrogen regulation vs UvrA
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

## Core review v0.3 (`core_review.py`, job 27330635, commit 06b11ea)

Scope: every `target_lt2500` core key with ≥ 1 different-context pair, its target-scope occurrences and its
different-context pairs. Cost < 1.5 min per subset. Evidence stays `assembly_only`; nothing here is an insertion
call or a site group.

**v0.2 boundary counts are superseded.** In v0.2 (job 27330567, `6b8e005`, reported in `45fa49f`), a chained
alignment only had to overlap the core, and a side that the alignment did not reach was recorded as extension 0.
A 100 bp local hit inside a 1,000 bp core therefore gave 0 / 0 and "resolved" (review counterexample, now a test).
v0.3 (`06b11ea`):
- measures a side only when the alignment covers that core end in both copies; otherwise it is unknown, never 0;
- counts the partners that do not cover an end;
- resolves a side only when every different-context partner gives a known extension;
- gives an estimated span only when both sides are resolved.

Occurrences resolved on both sides went from 219 / 314 / 169 to **123 / 214 / 122** (single-gene) and from
42 / 92 / 20 to **16 / 40 / 3** (two-gene). The sample split, identities and host context are unchanged.

What the review adds, per different-context pair and occurrence:
- **Sample split.** Pairs within one assembly (multi-copy) and across assemblies are both kept and counted
  separately. Within-assembly copies fit "same core, different context" but say nothing about strain variation.
- **Real similarity.** A family (MMseqs 0.9/0.8) is not an identical sequence. Each pair carries exact protein
  identity (sha256), blastp identity / coverage per core gene in core order (minimum for two-gene cores), and blastn
  identity / coverage of the oriented coding span.
- **Shared sequence beyond the CDS.** Each pair's coding span ± 3 kb is aligned (blastn, plus strand), and HSPs on
  the core diagonal are chained (gaps ≤ 50 bp). A side is measured only when the chain covers that core end in both
  copies. It is resolved when every partner is measured, none reaches the window or contig end, and the values
  differ by ≤ 20 bp. A resolved side is a **candidate conserved end**: sequence shared by copies in different
  contexts. It is not a complete insertion boundary. End hints (20 bp inverted-repeat identity, exact direct repeat
  ≤ 12 bp outside) are hints only.
- **Mobile-like neighbours.** A neighbour is mobile-like by product keyword (transposase / IS, integrase /
  recombinase / resolvase, phage) or because its family is itself a reviewed single-gene core. Host context = first
  K=3 neighbours per side after skipping mobile-like genes (scan ≤ 10 genes); pairs are re-related on it with the
  unchanged rule.
- **Genome ANI.** Cross-assembly pairs carry the skani ANI of the two assemblies (frozen S2 edges) and core DNA
  identity minus ANI. This is a diagnostic with no cutoff.

Single-gene cores (target scope):

| | comamonas | delftia | paracidovorax |
|---|---:|---:|---:|
| cores reviewed | 51 | 79 | 21 |
| with a cross-assembly different pair | 44 | 67 | 17 |
| within-assembly different pairs only | 7 | 12 | 4 |
| different pairs within / across assemblies | 188 / 411 | 215 / 541 | 669 / 1,287 |
| pairs with identical protein / identity < 95% | 93 / 135 | 272 / 202 | 790 / 130 |
| pairs with identical core DNA | 58 | 167 | 710 |
| occurrences | 242 | 374 | 200 |
| up side: resolved / no partner covers the end / some partners do not / inconsistent / reaches window | 203 / 15 / 10 / 14 / 0 | 270 / 56 / 13 / 35 / 0 | 151 / 27 / 4 / 18 / 0 |
| down side: same order | 142 / 77 / 9 / 14 / 0 | 257 / 69 / 22 / 23 / 3 | 154 / 25 / 5 / 16 / 0 |
| **both sides resolved** (mobile-keyword cores / other) | **123** (58 / 65) | **214** (49 / 165) | **122** (113 / 9) |
| resolved: estimated span < 2500 / ≥ 2500 | 123 / 0 | 214 / 0 | 122 / 0 |
| resolved: total shared extension (bp) | 1–635 | 2–1,804 | 12–272 |
| occurrences with a mobile-like near neighbour | 36 | 94 | 63 |

Two-gene cores (target scope):

| | comamonas | delftia | paracidovorax |
|---|---:|---:|---:|
| cores reviewed | 10 | 27 | 6 |
| with a cross-assembly different pair / within only | 8 / 2 | 26 / 1 | 6 / 0 |
| different pairs within / across assemblies | 14 / 44 | 7 / 92 | 2 / 15 |
| pairs with identical protein / identity < 95% | 4 / 22 | 2 / 50 | 0 / 5 |
| occurrences / both sides resolved | 42 / 16 | 97 / 40 | 20 / 3 |
| resolved: estimated span < 2500 / ≥ 2500 | 13 / 3 | 40 / 0 | 3 / 0 |

How to read it:
- **A resolved side is a candidate conserved end, not an insertion boundary.** Most unresolved sides have no
  partner alignment that covers the end. Their span is unknown, and neither a zero nor a fixed extension should be
  read into it.
- **IS families give consistent candidate ends.** For example, IS5 PF_4b70bb8e2b9d (Delftia) shares +62–64 / +9–10 bp
  beyond the CDS with every partner, and ISAav5 (Paracidovorax) shares +90–91 / +0 bp (end covered, measured 0). That
  fits the non-coding end of an element, but the boundary itself is not established. The IR hint is weak because
  BLAST ends drift a few bp, and 2–4 bp direct repeats are expected by chance.
- **RdgC and the ETF FixA/FixB cores are stable multi-copy paralogs, not single-copy genes.** The v0.2 text said
  "one copy per assembly"; that was wrong. RdgC (Comamonas) has 2–4 copies in each of 4 assemblies. The ETF case
  cores have exactly 2 copies in every assembly of each subset. Copies share 85–92% DNA (RdgC), and most of their
  ends are not covered by partner alignments (v0.3: unknown). "Same family, different neighbours" therefore also
  picks up paralogs that are present in every assembly at a constant copy number.
- **Other non-mobile cores** with one copy per assembly and 85–92% DNA identity (for example H-NS, LysR, permeases
  in Delftia) fit a gene in a rearranged neighbourhood between medium-distance assemblies, but this is not tested.
- **ANI comparison.** For assemblies at ANI ≥ 99 every copy is about 100% identical, so identity minus ANI says
  nothing. For more distant assemblies, conserved genes normally exceed the genome average. A per-pair background
  baseline can help **rank** unusually conserved cores. It cannot by itself tell movement from paralogy or from
  inheritance followed by rearrangement.
- **No nested case was identified by the current rule.** That does not mean there is no nesting. Mobile-like near
  neighbours are common (36 / 94 / 63 single-gene occurrences), but no different-context pair became host-similar
  after skipping them. Pairs only moved to insufficient labels or one side shared (single-gene 3 / 8 / 57 pairs;
  two-gene 0 / 1 / 0). The "reviewed core family" criterion also flags ordinary paralogs (for example SGNH
  hydrolase), with 5 / 42 / 16 occurrences (single- and two-gene) flagged by it alone, so it is a broad flag.

Representative cases (fixed rule in `review_summary.json`: A = most cross-assembly host-different pairs, B = most
within-assembly pairs with none across, D = two-gene; no case C because the current rule found none; figure
`core_review/cases_review.svg`, rows with the partner counts in `core_review/case_occurrences.tsv`):

| subset | case | core | occ / samples | different pairs within / across | protein id % | core DNA id % | shared beyond CDS, drawn copies (up / down) |
|---|---|---|---|---|---|---|---|
| comamonas | A1 | PF_147add17340c ISL3 | 29 / 3 | 102 / 221 | 98.3–100 | 98.9–100 | +83–148 / +3–4; one down side inconsistent |
| comamonas | A2 | PF_e2e705caf61f RdgC (2–4 copies per assembly) | 12 / 4 | 8 / 24 | 89.6–91.9 | 85.0–89.3 | +1 / unknown (no partner covers the end) |
| comamonas | B | PF_8e9360ea03b6 IS5 | 5 / 1 | 10 / 0 | 100 | 99.9–100 | +57–58 / +18 |
| comamonas | D | ETF alpha/FixB + partner (2 per assembly) | 10 / 5 | 5 / 20 | 90.9–99.4 | 86.5–99.4 | unknown / not covered by some partners |
| delftia | A1 | PF_4b70bb8e2b9d IS5 | 30 / 4 | 44 / 134 | 99.3–100 | 97.8–100 | +62–64 / +9–10 |
| delftia | A2 | PF_3070024975f6 IS481 | 21 / 5 | 45 / 112 | 98.2–100 | 98.4–100 | +171–206 / +227–233, mostly inconsistent |
| delftia | B | PF_ad056eba1256 IS5 | 15 / 1 | 44 / 0 | 100 | 100 | +46–48 / +9–12 |
| delftia | D | ETF alpha/FixB + partner (2 per assembly) | 14 / 7 | 6 / 30 | 94.2–99.6 | 87.6–98.8 | unknown / unknown |
| paracidovorax | A1 | PF_34809b3c6b4f ISAav5 | 54 / 3 | 350 / 701 | 99.4–100 | 99.7–100 | +90–91 / +0–1 |
| paracidovorax | A2 | PF_f094f2564f99 IS5 | 30 / 4 | 75 / 212 | 92.2–100 | 86.7–100 | +56–57 / +9–10 |
| paracidovorax | B | PF_282ac9979b1c IS5 | 8 / 2 | 20 / 0 | 100 | 100 | +67–68 / +133–134 |
| paracidovorax | D | ETF alpha/FixB + partner (2 per assembly) | 8 / 4 | 2 / 6 | 95.8–98.4 | 94.1–98.7 | not covered by some partners on both sides |

Identity ranges are over the case's different-context pairs; extensions are over the drawn occurrences.

Outputs per subset in `<subset>/protein_entry/core_review/`: `review_summary.json`, `review_cores.tsv`,
`review_occurrences.tsv.gz` (with partners not covering each end), `review_pairs.tsv.gz`, `case_occurrences.tsv`,
`cases_review.svg`.

Still open:
- complete core boundaries: the unknown sides, IR / target-site checks on aligned ends rather than fixed 20 bp
  windows, and graph-block merging only under consistent adjacency, order and orientation;
- a per-pair background baseline (similar-context single-copy genes between the same two assemblies), for ranking
  only;
- separating paralogs present in every assembly at a constant copy number from copies whose number varies;
- a narrower mobile-neighbour rule (the core-family criterion over-flags paralogs);
- site groups (S6) and read checks.

Evidence stays `assembly_only`, `production_ready=false`.
