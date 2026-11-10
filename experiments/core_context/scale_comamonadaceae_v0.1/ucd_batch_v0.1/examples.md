# Real U–X–D examples (`candidate_insertion`, final_status = pass)

One per pair scope × structure, chosen by a fixed rule: the longest insertion, then the smallest pair_id (`ucd_batch_examples.py`).
- Sequences come from the Check-1 evidence: re-extracted from the verified source FASTA and identical to the screen's, base by base.
- HSPs are those of Check 2. `chain` marks the collinear HSPs that make up the insertion match.
- Coordinates are 0-based half-open on the contig forward strand. Sequences are shown in the insertion orientation. B is shown as read; when `b_flipped` is true, B's insertion is the reverse complement of A's.
- This is a structural pattern only: it says nothing about mobility or an empty site.

## cross_assembly / single_gene: `M0507f5b9f4d3`

core `PF_4090266a9173`, core match `identical_core_dna` (core DNA identity 100.0), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `s|SAMEA5159697|NZ_LR594660.1|560` | GCF_901827155.1 | NZ_LR594660.1 (839194) | + | 533153-534153 | 534153-536654 | 536654-537654 | `49c6d52846acbea3…` |
| B | `s|SAMN31283087|NZ_CP109857.1|3202` | GCF_036828685.1 | NZ_CP109857.1 (6759232) | + | 3517648-3518648 | 3518648-3521149 | 3521149-3522149 | `741f8ddb8982515f…` |

- **Insertion:** X/X′ identity 99.84%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2501 / 2501 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  cggggcacctggcctacttcatcagcgacacggtcgatgc | GCTGTAAGCGTTCCGCGAACCCATCTTTCATCTCACCGTT ... GGCGCAAGGGTGTCGTTCGCGGAACGCTTACGGTGTGCCG | gcggcgccggcgcagttgctgccggtgcgaatcgagccga
B  ccggctctttgggccggttttttttcgcctgcggttttct | GCCGTAAGCGTTCCGCGAACCCATCTTTCATCTCACCGTT ... GGCGCAAGGGTGTCGTTCGCGGAACGCTTACTTTCTGCCG | atgtcggcatgaaaaaaccggccgcccttgcgggccggcc
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 99.84 | 2501 | 0–2501 | 0–2501 | yes |
| X | X | minus | 94.118 | 34 | 819–853 | 819–853 |  |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

## cross_assembly / two_gene: `M6049f5a0f6cf`

core `PF_4090266a9173+,PF_e82e9faa7ca9+`, core match `identical_core_dna` (core DNA identity 100.0), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `t|SAMEA5159699|NZ_LR594667.1|340` | GCF_901827165.1 | NZ_LR594667.1 (801406) | - | 301415-302415 | 298914-301415 | 297914-298914 | `8da96a181e4970e8…` |
| B | `t|SAMN31283087|NZ_CP109857.1|3202` | GCF_036828685.1 | NZ_CP109857.1 (6759232) | + | 3517648-3518648 | 3518648-3521149 | 3521149-3522149 | `741f8ddb8982515f…` |

- **Insertion:** X/X′ identity 99.84%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2501 / 2501 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  cggggcacctggcctacttcatcagcgacacggtcgatgc | GCTGTAAGCGTTCCGCGAACCCATCTTTCATCTCACCGTT ... GGCGCAAGGGTGTCGTTCGCGGAACGCTTACGGTGTGCCG | gcggcgccggcgcagttgctgccggtgcgaatcgagccga
B  ccggctctttgggccggttttttttcgcctgcggttttct | GCCGTAAGCGTTCCGCGAACCCATCTTTCATCTCACCGTT ... GGCGCAAGGGTGTCGTTCGCGGAACGCTTACTTTCTGCCG | atgtcggcatgaaaaaaccggccgcccttgcgggccggcc
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 99.84 | 2501 | 0–2501 | 0–2501 | yes |
| X | X | minus | 94.118 | 34 | 819–853 | 819–853 |  |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

## same_assembly / single_gene: `Md8105b6b875b`

core `PF_003edf97e510`, core match `identical_core_dna` (core DNA identity 100.0), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `s|SAMN57441782|NZ_CM170111.1|3133` | GCF_057218865.1 | NZ_CM170111.1 (5140579) | + | 3331554-3332554 | 3332554-3335062 | 3335062-3336062 | `72d6ebbbdc043976…` |
| B | `s|SAMN57441782|NZ_CM170111.1|348` | GCF_057218865.1 | NZ_CM170111.1 (5140579) | - | 379482-380482 | 376971-379482 | 375971-376971 | `72d6ebbbdc043976…` |

- **Insertion:** X/X′ identity 99.8%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2508 / 2511 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  gttgttgggtcgatgcgatacagcgtggacctggcctgac | TGTTGATTTCCACCCAGAAGTGACCCACTAGCCCCCGGTT ... GCGGGTCAGTTCTGCGTGGAAATCAACACAAAGTAGGCAA | gaagctacgcccgctttgggcaggcccacattgccagagc
B  tggcttcgctgcacagcctttccctaccaaggccttgaag | TGTTGATTTCCACCCAGAAGTGACCCACTAGCCCCCGGTT ... GGTCAGTTCTGCGTGGAAATCAACACCTTGAAGGAGGCAA | cggttggacgcagtgccaataacgcaccatttctgcttgc
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 99.801 | 2511 | 0–2508 | 0–2511 | yes |
| X | X | minus | 75.806 | 62 | 2444–2496 | 0–62 |  |
| X | X | minus | 75.806 | 62 | 0–62 | 2444–2496 |  |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

## same_assembly / two_gene: `Mf47160818cc6`

core `PF_003edf97e510+,PF_7730b136b319+`, core match `identical_core_dna` (core DNA identity 100.0), s6_grouping label `grouping_unresolved` (not used)

| side | occurrence | assembly | contig (length) | strand | U | X | D | FASTA sha256 |
|---|---|---|---|---|---|---|---|---|
| A | `t|SAMN57441782|NZ_CM170111.1|3133` | GCF_057218865.1 | NZ_CM170111.1 (5140579) | + | 3331554-3332554 | 3332554-3335062 | 3335062-3336062 | `72d6ebbbdc043976…` |
| B | `t|SAMN57441782|NZ_CM170111.1|347` | GCF_057218865.1 | NZ_CM170111.1 (5140579) | - | 379482-380482 | 376971-379482 | 375971-376971 | `72d6ebbbdc043976…` |

- **Insertion:** X/X′ identity 99.8%, coverage 1.0 / 1.0 (collinear chain of 1 HSP, plus strand); lengths 2508 / 2511 bp.
- **Contexts** (coverage A / B by HSPs ≥ 80%): UU different (0.0 / 0.0); DD different (0.0 / 0.0); UU_inv different (0.0 / 0.0); DD_inv different (0.0 / 0.0); UD different (0.0 / 0.0); UD_inv different (0.0 / 0.0); DU different (0.0 / 0.0); DU_inv different (0.0 / 0.0)
- **Verdicts:** screen `ucd_supported`; source_check `pass`; ucd_check `ucd_supported`.

Junctions (40 bp of U | 40 bp of X … 40 bp of X | 40 bp of D):

```
A  gttgttgggtcgatgcgatacagcgtggacctggcctgac | TGTTGATTTCCACCCAGAAGTGACCCACTAGCCCCCGGTT ... GCGGGTCAGTTCTGCGTGGAAATCAACACAAAGTAGGCAA | gaagctacgcccgctttgggcaggcccacattgccagagc
B  tggcttcgctgcacagcctttccctaccaaggccttgaag | TGTTGATTTCCACCCAGAAGTGACCCACTAGCCCCCGGTT ... GGTCAGTTCTGCGTGGAAATCAACACCTTGAAGGAGGCAA | cggttggacgcagtgccaataacgcaccatttctgcttgc
```

Check-2 HSPs (q = A, s = B as read):

| q | s | strand | identity | length | q range | s range | chain |
|---|---|---|---|---|---|---|---|
| X | X | plus | 99.801 | 2511 | 0–2508 | 0–2511 | yes |
| X | X | minus | 75.806 | 62 | 2444–2496 | 0–62 |  |
| X | X | minus | 75.806 | 62 | 0–62 | 2444–2496 |  |

No HSP other than X/X′ means no reliable similarity between any pair of flanks.

