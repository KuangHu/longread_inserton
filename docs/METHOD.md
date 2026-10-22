# Fixed method and interpretation

## Native discovery

Each sample is an assembly set (chromosomes and plasmids of one isolate).
Samples are never pooled into one discovery search. Commands are unchanged from
the native-only pilot:

```bash
pals -self INPUT -out hits.gff -log pals.log
piler2 -trs hits.gff -out families.gff -piles piles.gff -images images.gff -log piler.log
piler2 -trs2fasta families.gff -seq INPUT -path fams
muscle -in fams/FAMILY -out aligned_fams/FAMILY -maxiters 1 -diags1
piler2 -cons aligned_fams/FAMILY -out cons/FAMILY -label FAMILY

build_lmer_table -sequence INPUT -freq lmers.freq
RepeatScout -sequence INPUT -freq lmers.freq -output families.fasta -ranges ranges.tsv
```

No BLAST-to-GFF replacement, filter-stage-1/2, custom copy/length/complexity
rejection, RepeatMasker, TIR or TSD gate is used. Native defaults have their own
constraints: PALS defaults to 400 bp / 94% identity and PILER-DF requires a
dispersed family. PILER's default `-trs2fasta` member cap is 16; it can affect
the members contributing to native consensus. The inventory and boundary
analysis use **all GFF records**, not just the extracted subset. A native family
without a consensus is retained and explicitly marked; BLASTX cannot annotate
that missing consensus. Native duplicate coordinates remain duplicate records.

References: [PILER manual](https://www.drive5.com/piler/piler_userguide.html),
[PALS manual](https://www.drive5.com/pals/pals_userguide.htm),
[maintained RepeatScout](https://github.com/Dfam-consortium/RepeatScout).

## Post-hoc annotation

All provided EMBL/GenBank records must match their assembly sequences exactly
(case-insensitive). CDS overlap uses every location part, never the outer
bounding interval of an origin-crossing gene. Protein translations come from
deposited annotations; no ORF predictor or HMM is substituted.

Family consensuses are searched with BLASTX against the pooled deposited
proteins of the samples supplied in this run, `E <= 1e-5`,
`max_target_seqs=1000`. All emitted HSPs remain in `deposited.blastx.tsv`.
The `tnp-supported` subset requires both:

1. at least one significant protein hit whose title contains `transposase`;
2. at least one original member overlapping a CDS whose product contains that word.

These criteria reproduce the current pilot, including its dependence on
annotation naming and the available protein set. They are not proof of active
transposition. Support is at family level; every member of a selected family is
retained, including members without a complete annotated tnp. No protein hit is
not evidence of absence. `--subset all` removes only this post-hoc selection;
it does not change discovery. The full Swiss-Prot comparison from the exploratory
work is not a requirement for this package and is not rerun here.

## Independent homology extent

For each selected valid native member:

1. Extract the native interval plus 2000 bp on each side, clipped at contig ends;
   orient according to the native family strand.
2. Compare all contexts using BLASTN (`task=blastn`, `word_size=11`, `E=1e-5`,
   `dust=yes`, `soft_masking=true`, `max_target_seqs=100000`, BTOP output).
3. Within each directed pair of different members in the same family, consider
   forward HSPs in that oriented frame. Count BTOP alignment columns where
   **both** positions are inside their respective native cores. Require at least
   one linked base; choose the HSP with most linked bases, then highest bitscore.
   This is an HSP association rule, not a candidate rejection threshold.
4. For each member, take the lower observed median of starts and, separately,
   ends over its peer HSPs. Save all peer observations, endpoint ranges and
   independent `inferred_start0/end0`.

Linking both cores avoids mistaking overlapping physical context windows for
two aligned repeat copies. The median is a provisional homology-block estimate,
not a statistical confidence interval. Low coverage, fragments, shared external
sequence and heterogeneous extensions can still make it an inaccurate full-element
boundary. The two medians can come from different peer comparisons.

## Reference annotation and export status

All deposited mobile_element annotations overlapping the native interval are
listed. Export adopts an annotation only if exactly one contiguous feature
contains the native interval midpoint. Origin-spanning compound mobile_element
features are reported but not represented by a chromosome-wide bounding box.
This association rule prevents a 2-bp overlap with a neighboring IS from adopting
that entire neighboring element. Independent inferred coordinates remain available.

| `boundary_source` | Exported range |
|---|---|
| `deposited_mobile_element` | Existing annotation; inference recorded separately if available |
| `provisional_homology_median` | Median peer homology extent; not a certified complete insertion |
| `homology_extent_censored_at_window` | Median extent reaches a context edge; incomplete extent |
| `native_interval_only_unresolved` | Original native range; no peer HSP linking both cores |
| `invalid_native_coordinates` | Inventory only; no fabricated sequence export |

Any individual peer HSP touching the context edge is additionally flagged even
if the median does not touch. A candidate with annotated boundaries may have no
independent peer support. A family with only invalid members remains in the
selected-family inventory, although it has no sequence-based summary.

## Flanks and comparisons

Export 500 bp per side when available. Negative-strand records are reverse
complemented for the family-oriented triplet; separate genomic-forward exports
are included. Flanks can be truncated or empty at contig ends, without deleting
the record or wrapping the origin. Empty flank entries are kept in the exports
but omitted from the BLAST query FASTA.

Same-family up-up and down-down flanks are compared with the same BLASTN settings.
Per-pair observations include the HSP with longest minimum span on the two
sequences and any physical overlap of the genomic flank windows. The 250-bp
summary is descriptive, not a filtering rule. A missing HSP is uninformative
when a flank is absent and is never proof of unrelated ancestry. No up-down
rearrangement interpretation is attempted.

Cross-tool overlaps are all pairs of selected native PILER/RS records sharing
bases on the same sample and contig. Exported boundary differences are reported;
there is no automatic locus deduplication or assertion of equivalent full elements.
Fixed-boundary direct-repeat (2–20 bp) and terminal inverted-repeat (4–50 bp)
lengths are descriptive strings, not validated TSD/TIRs or boundary search criteria.

The native discovery engine already favors certain repeat geometries, especially
PILER-DF isolated copies. External-flank divergence is therefore supporting
description, not a fully independent validation of discovery specificity.
