# Real U–X–D examples (`candidate_insertion`, final_status = pass)

One per pair scope × structure, chosen by a fixed rule: the longest insertion, then the smallest pair_id (`ucd_batch_examples.py`).
- Sequences come from the Check-1 evidence: re-extracted from the verified source FASTA and identical to the screen's, base by base.
- HSPs are those of Check 2. `chain` marks the collinear HSPs that make up the insertion match.
- Coordinates are 0-based half-open on the contig forward strand. Sequences are shown in the insertion orientation. B is shown as read; when `b_flipped` is true, B's insertion is the reverse complement of A's.
- This is a structural pattern only: it says nothing about mobility or an empty site.

## cross_assembly / single_gene: `M4c2b936d3f6b`

core `PF_579bd8cde7ea`, core match `homolog_full_coverage` (core DNA identity 95.128), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `s|SAMN02472164|NZ_AFPT01000135.1|104` | GCF_000218805.1 | NZ_AFPT01000135.1 (131518) | - | 112810-113810 | 110311-112810 | 109311-110311 | `d676e45e98c01b01…` |
| B | `s|SAMN05714654|NZ_CP028287.1|3826` | GCF_003029805.1 | NZ_CP028287.1 (5574669) | - | 4870832-4871832 | 4868333-4870832 | 4867333-4868333 | `48b31c47b2f49582…` |

- **Insertion:** X/X′ identity 95.28%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2499 / 2499 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  ccgtcgtgctggaaggggtggccggggcgtccttgccggg | GCCGCTCGCCGTCGCCGAGCCGGGCGGGGCGCCGGTGGCC ... GCGCCAGCCAGCCCGATCCCGCTCCCCGTCCCGCACCGCC | ccgcgcgggcggtccccggcccgcgggccgcttcaccctg
B  ggcacccggggctgccggtgggggctcccgccctggccca | GCCGCTCGCCGTCGCCGAGCCGGGTGGGGCGCCGGTGGCC ... GCGCCAGCCAGCCCGATCCCGCTCCCCGTTCCGCACCGCC | ggtgagccggtcggtgtcggnnnnnnnnnnnnnnnnnnnn
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 95.278 | 2499 | 0–2499 | 0–2499 | yes |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

## cross_assembly / two_gene: `Ma94c20621031`

core `PF_0b409518c7ff+,PF_656b7e14456e+`, core match `homolog_full_coverage` (core DNA identity 99.957), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `t|SAMN43282391|NZ_JBHUMV010000006.1|93` | GCF_042685145.1 | NZ_JBHUMV010000006.1 (390816) | - | 103471-104471 | 100973-103471 | 99973-100973 | `170a9a4ba058c5a1…` |
| B | `t|SAMN56682558|NZ_CM151366.1|74` | GCF_056488365.1 | NZ_CM151366.1 (3821750) | + | 85865-86865 | 86865-89363 | 89363-90363 | `8678b8195daeed9a…` |

- **Insertion:** X/X′ identity 99.92%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2498 / 2498 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  ggcaacaagggcgagcgcgcaatccacttttttctgcggc | GAGTGTTGATTTCCACGCAGAAGTGACCCACTAACCCCCA ... GAGATGGAATTGGTGGGTTCAGATTGCATGGAAATCAACA | gcggcgagctccacggcttccttgaggatctcgttctcct
B  ttagccaaaaccgagtcctgatctcatcaatctgttcctg | GATTGTTGATTTCCACGCAGAAGTGACCCACTAACCCCCA ... GAGATGGAATTGGTGGGTTCAGATTGCATGGAAATCAACA | ttcctggatagacatctgcaaaagcttccaaagcgcaaaa
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 99.92 | 2498 | 0–2498 | 0–2498 | yes |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

## same_assembly / single_gene: `M0261582d02af`

core `PF_003edf97e510`, core match `identical_core_dna` (core DNA identity 100.0), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `s|SAMN28016767|NZ_CP096918.1|1536` | GCF_027595065.1 | NZ_CP096918.1 (4015322) | + | 1586697-1587697 | 1587697-1590196 | 1590196-1591196 | `d521108554b39cb0…` |
| B | `s|SAMN28016767|NZ_CP096918.1|8` | GCF_027595065.1 | NZ_CP096918.1 (4015322) | - | 12241-13241 | 9742-12241 | 8742-9742 | `d521108554b39cb0…` |

- **Insertion:** X/X′ identity 100.0%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2499 / 2499 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  aattaccttgaccagctcaagcagggcgccgaggtcattc | AGTGTTGATTTCCATCCAGAAGTGACCCACTAGCCCCCGG ... AGATGGAAATTGCGGGTCAGTTCTGCGTGGAAATCAACAC | actggatacttattcccgtaggacgcatagatggatgtac
B  gccactcagcttggcatgttagcgctctctttgcatgcgg | AGTGTTGATTTCCATCCAGAAGTGACCCACTAGCCCCCGG ... AGATGGAAATTGCGGGTCAGTTCTGCGTGGAAATCAACAC | tcgcttctccgggtgttcctgctttgtcgaggtttgcgca
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 100.0 | 2499 | 0–2499 | 0–2499 | yes |
| X | X | minus | 74.603 | 63 | 2446–2499 | 1–64 |  |
| X | X | minus | 74.603 | 63 | 1–64 | 2446–2499 |  |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

## same_assembly / two_gene: `M1deff05b6c1c`

core `PF_003edf97e510+,PF_7730b136b319+`, core match `homolog_full_coverage` (core DNA identity 99.913), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `t|SAMN07626132|NZ_PEFL01000001.1|754` | GCF_002754475.1 | NZ_PEFL01000001.1 (3129460) | + | 779631-780631 | 780631-783130 | 783130-784130 | `cc7b9947e7d67c98…` |
| B | `t|SAMN07626132|NZ_PEFL01000003.1|575` | GCF_002754475.1 | NZ_PEFL01000003.1 (774302) | - | 625257-626257 | 622758-625257 | 621758-622758 | `cc7b9947e7d67c98…` |

- **Insertion:** X/X′ identity 99.88%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2499 / 2499 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  cagcaatccatctgcctgcacttgttgcccactgtcagct | GCTTGTTGATTTCCACCCAGAAGTGACCCACTAGCCCCCG ... GAGATGGAAATTGCGGGTCAGATTTGCGTGGAAATTAACA | gggtgtcgatcacaacaagggcgggtggtggtgtgtattc
B  ttgcttgtggtcatgggttgtattgataaccggtcactca | GCGTGTTGATTTCCACCCAGAAGTGACCCACTAGCCCCCG ... GAGATGGAAATTGCGGGTCAGATTTGCGTGGAAATTAACA | ctgagtcactcatacgatttttcatttattgtttcatggc
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 99.88 | 2499 | 0–2499 | 0–2499 | yes |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

