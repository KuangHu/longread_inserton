# S6 site groups v0.1 (prototype rules)

Code: `site_groups.py` (tests: `test_site_groups.py`). Inputs: one subset's committed `protein_entry/` and
`protein_entry/core_review/`, plus the scratch protein work dir (genes, MMseqs assignments). Prototype on the three
frozen subsets only (`subsets_20261027`). Evidence stays `assembly_only`; `production_ready=false`.

## What a site group is

A site group is **repeated observation of one core architecture at a corresponding genome position**. It is not an
independent insertion event. Gap lists and empty sites are neither used nor required.

```
similar core + corresponding left and right host context  -> same site group (candidate)
similar core + evidence for a different position          -> different site groups (candidate)
insufficient context or unclear correspondence            -> unresolved
```

- **Architecture.** Single-gene and two-gene cores, and different gene compositions, are never mixed: grouping is
  per `structure + core_key`.
- **Occurrences.** Every `target_lt2500` occurrence of every core reviewed by `core_review.py` (cores with ≥ 1
  different-context pair). This includes occurrences without a different-context partner.

## Context windows

Each window is per side, in core orientation:
- **near:** the first K = 3 genes as they are, mobile or not;
- **host:** the first N = 5 genes after skipping genes whose product matches a mobile keyword (transposase / IS,
  integrase / recombinase / resolvase, phage). At most 15 genes are scanned. A family being a candidate core is never
  a reason to skip.

Labels are homology groups (CDS with protein) or RNA type + product. CDS without a protein have no label. Both
windows are recorded per occurrence, along with the number of skipped genes.

## Side relation (host windows of two occurrences)

- Labels occurring more than once in either window are ambiguous and never count as matches (repeated flanks).
- A match needs the same label and the same direction relative to the core.
- Matched labels must keep the same order in both windows; the support is the longest increasing chain.
- `same`: support ≥ 2.
- `different`: both windows full (5 genes, all labelled) and no shared label.
- `unknown`: everything else. Missing labels, short windows (contig ends), or shared labels out of order or
  orientation are never `different`.

The old one-shared-label-in-three-neighbours relation is not used for grouping.

## Pair verdict

| verdict | rule |
|---|---|
| `different_position_same_assembly` | same assembly, non-overlapping; if both sides are also `same`, the pair carries `repeated_flank` |
| `same_position_candidate` | different assemblies, both sides `same` |
| `different_position_candidate` | different assemblies, both sides `different` |
| `conflicting_sides` | one side `same`, the other `different` |
| `unresolved` | anything else (including overlapping occurrences in one assembly) |

**Multiple correspondence.** If an occurrence has `same_position_candidate` links to two or more occurrences of one
other assembly, those links are kept as `ambiguous_multiple_correspondence` and are not used for grouping.

## Grouping (no transitive merging)

1. Take the remaining `same_position_candidate` links and form connected components only as *candidates*.
2. Accept a component as a site group only if **every** pair inside it is `same_position_candidate` (a clique) and
   **no assembly appears twice**. `A≈B, B≈C` with `A` vs `C` not same, or two distinct positions of one assembly
   reached through other samples, rejects the component. All its members stay `unresolved` with reason
   `conflict_in_component`.
3. An occurrence with no same-position link whose every pair is a different-position verdict is a
   `single_observation_site`. That covers a position seen once, with evidence that it differs from every other
   observation of the core.
4. Everything else is `unresolved`, with reason `ambiguous_multiple_correspondence`, `conflicting_sides` or
   `unresolved_pairs`. **Unresolved occurrences are never counted as sites.**

Every accepted group lists its members, locations and per-pair evidence (`s6_group_evidence.tsv`: ordered matches per
side, repeated-label flag).

## Stable-end list

Built from `core_review` outputs, for occurrences with ≥ 1 different-context partner:
- **Spans:** coding span, shared-sequence estimate (coding span + median shared extension) and per-side status.
- **Partners:** valid partners per side and on both sides, the assemblies they come from, the site groups they
  belong to (unresolved partners counted separately), and a `single_partner` flag.
- **Length:** a strict `< 2500 bp` re-check on the shared estimate when both sides are resolved.
- **Two lists:** `stable_end_candidates.tsv` (both sides resolved: candidate conserved ends, not confirmed insertion
  boundaries) and `boundary_review_unknown.tsv` (any unknown side). No candidate is dropped.

## Reported, never used to exclude

Copy number per assembly (0 = no instance meeting the current complete-CDS criteria was detected; not biological
absence), product keywords, near mobile neighbours, protein / DNA identity and genome ANI are reporting and ranking
fields only.

## Outputs (`<subset>/protein_entry/s6/`)

- `s6_occurrences.tsv.gz`: every occurrence with coordinates, strand, coding span, boundary status, near and host
  windows, skipped mobile genes, S6 status, site group or unresolved reason.
- `s6_pairs.tsv.gz`: every pair within a core, with side relations, supports, verdict and flags.
- `s6_site_groups.tsv`, `s6_group_evidence.tsv`, `s6_cores.tsv`.
- `stable_end_candidates.tsv`, `boundary_review_unknown.tsv`.
- `s6_summary.json`, `cases_s6.svg`.

## Known limits of v0.1

- Gene-order evidence only. No flank DNA alignment and no empty-site evidence.
- A strict clique rule: one unresolved pair blocks a whole component. Some real groups are left unresolved by design.
- Different-position evidence between distant assemblies can also come from rearrangement around a conserved site.
- Mobile skipping is by product keyword, so unannotated mobile DNA stays in the host window.

## Review layer (`s6_review.py`, tests `test_s6_review.py`)

This layer reads the outputs above and never changes them. It writes to `<subset>/protein_entry/s6/review/`.

**Boundary denominators.** Every occurrence is in exactly one class:
- `both_sides_resolved`;
- `reviewed_side_unknown`;
- `not_reviewed_no_different_partner`: core_review never measured its ends.

Stable ends are split into:
- one valid partner;
- ≥ 2 partners in one assembly;
- partners from ≥ 2 assemblies;
- partners in ≥ 2 position units.

One valid partner means one comparison supports that end.

**Coordinates vs host context.** A non-overlapping pair in one assembly is always at different coordinates. Its
host context is classified separately from the side relations:

| class | sides |
|---|---|
| `context_different_both_sides` | both `different`: the target U–core–D pattern |
| `context_repeated_both_sides` | both `same` (`repeated_flank`) |
| `context_repeated_one_side` | one `same`, one `different` |
| `context_unknown` | any side `unknown` |

Flank DNA cross-check, reported only and never used to change a class:
- blastn between the two copies' outer flanks, 500–3,000 bp beyond each core end in core orientation;
- the inner 500 bp is left out because element ends can extend past the CDS;
- coverage is the smaller of the two copies' covered fractions;
- classes: `shared` ≥ 0.5, `not_shared` ≤ 0.1, otherwise `partial`, and `short_flank` if the outer flank is < 1,000 bp.

**Position units and their relations.** A unit is a site group or a single-observation site. Unresolved occurrences
are not units. Two units of one core are compared over every cross pair of their members:
- `ambiguous`: some cross pair is `same_position_candidate`. This can only be a multiple-correspondence link.
- `conflicting`: some cross pair has `conflicting_sides`.
- `different`: every cross pair is a different-position verdict. The relation also carries a context: context
  different on both sides in every pair, repeated context in some pair, or unknown context in some pair.
- `unknown`: anything else. In particular, two accepted groups are **not** different by default.

Per core the report gives:
- the largest set of units that are pairwise `different`;
- the same with context different on both sides (exact maximum clique);
- never the plain unit count as a number of positions.

**Representative cases.** Fixed rules, single-gene only:
- context different: a same-assembly pair of supported positions;
- repeated context: a same-assembly pair with `repeated_flank`;
- conserved paralog positions: two non-mobile site groups, related `different` with context different, that share
  assemblies;
- multiple correspondence.

If no core meets a rule in a subset, that case is reported as missing.

## Flank DNA re-check v2 (`flank_dna.py`, tests `test_flank_dna.py`)

It replaces the v0.1 DNA cross-check of the review layer. It is auxiliary only.
- Both copies' outer windows (500–3,000 bp beyond each core end, in core orientation) are cut **before** alignment.
- Same side forward, same side reverse and cross side on each strand are recorded separately. The class comes from
  same side forward.
- It does not use or change the shared `core_review.blast`.
- `copy_gap_bp`, `windows_can_overlap` (gap < 6 kb) and `host_genes_shared_physically` flag neighbouring copies whose
  windows are the same DNA.
- `not_shared` = no outer-flank alignment detected by this check.
