# Pool refresh 2026-11-08

Run on lrc-xfer by `refresh_20261108.sh` (build_pool fetch + build assembly_first_v1, then download_fna on
the existing family data roots: only unverified rows are downloaded). `diff/` is `pool_refresh.py diff` against
the frozen pools (`scaleup_stage1_20261026`).

| family | frozen pool | refreshed pool | unchanged | new BioSample | updated accession | removed | new long-read runs |
|---|---|---|---|---|---|---|---|
| microbacteriaceae | 3444 | 3453 | 3443 | 9 | 1 | 0 | 2 |
| actinomycetaceae | 1139 | 1140 | 1139 | 1 | 0 | 0 | 0 |
| flavobacteriaceae | 3704 | 3717 | 3695 | 13 | 9 | 0 | 1 |
| eggerthellaceae | 329 | 329 | 328 | 0 | 1 | 0 | 0 |
| lachnospiraceae | 3254 | 3256 | 3254 | 2 | 0 | 0 | 0 |
| comamonadaceae | 1705 | 1716 | 1701 | 11 | 4 | 0 | 2 |
| nitrobacteraceae | 1975 | 1980 | 1975 | 5 | 0 | 0 | 0 |
| micrococcaceae | 2975 | 2995 | 2971 | 20 | 4 | 0 | 3 |
| metamycoplasmataceae | 1738 | 1760 | 1738 | 22 | 0 | 0 | 1 |
| erwiniaceae | 2500 | 2506 | 2495 | 6 | 5 | 0 | 0 |
| total | 22763 | 22852 | 22739 | 89 | 24 | 0 | 9 |

- All gap downloads verified. Failed rows in the s1 reports (Lachnospiraceae GCF_902652075.1, GCF_902652265.1;
  Erwiniaceae GCA_003568535.1) failed the same way in the frozen pool and are not in any catalogue manifest.
- U-core-D catalogue: 67 of the 89 new BioSamples are eligible (22 lack GFF3/FAA); they were added to batch 8 of
  `core_context/scale_joint_v0.2/manifest_b1-8r.tsv` (see `manifest_b1-8r.json`). The commit message of d6e781b
  says 92 new BioSamples and 25 without GFF3/FAA; the correct numbers are 89 and 22.
- The 24 updated accessions keep their frozen version in this round; replacing them would change pairs of
  earlier batches and is left for a later round.
