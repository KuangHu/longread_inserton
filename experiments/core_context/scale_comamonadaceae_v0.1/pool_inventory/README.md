# Wave-1 FNA pools: old pipeline vs the current U–core–D catalogue (`../../pool_inventory.py`)

`assemblies.tsv.gz` has one row per pool assembly (BioSample + accession.version) of the 10 wave-1 families. It
records:
- the old FNA / gap pipeline stages: s1 download and verify, s2 joint FNA and skani, s3 gap-caller pairs, s4 locus
  annotation;
- whether the assembly is in the current U–core–D catalogue;
- how many of its pair rows ran the two current checks (`ucd_batch_v0.2`), passed, and are in the strict
  candidate table.

| family | pool assemblies | s1 verified | s2 | s3 gap pairs | s4 loci | current U-core-D status |
|---|---:|---:|---:|---:|---:|---|
| Actinomycetaceae | 1,139 | 1,139 | 1,139 | 1,116 | 1,097 | old_pipeline_only 1,139 |
| Comamonadaceae | 1,705 | 1,705 | 1,705 | 1,691 | 1,617 | ucd_checks_done 1,679, ucd_pool_not_included 23, ucd_catalogue_no_pair 3 |
| Eggerthellaceae | 329 | 329 | 329 | 325 | 312 | old_pipeline_only 329 |
| Erwiniaceae | 2,500 | 2,499 | 2,499 | 2,486 | 2,361 | old_pipeline_only 2,500 |
| Flavobacteriaceae | 3,704 | 3,704 | 3,704 | 3,539 | 3,266 | old_pipeline_only 3,704 |
| Lachnospiraceae | 3,254 | 3,252 | 3,252 | 3,111 | 3,035 | old_pipeline_only 3,254 |
| Metamycoplasmataceae | 1,738 | 1,738 | 1,738 | 1,721 | 1,410 | old_pipeline_only 1,738 |
| Microbacteriaceae | 3,444 | 3,444 | 3,444 | 3,424 | 3,195 | old_pipeline_only 3,444 |
| Micrococcaceae | 2,975 | 2,974 | 2,974 | 2,963 | 2,889 | old_pipeline_only 2,975 |
| Nitrobacteraceae | 1,975 | 1,975 | 1,975 | 1,971 | 1,958 | old_pipeline_only 1,975 |

- **Status values:**
  - `ucd_checks_done`: a side of at least one catalogue pair whose two checks ran.
  - `ucd_catalogue_no_pair`: in the catalogue manifest, but no core of the assembly pairs with another occurrence.
  - `ucd_pool_not_included`: left out of the catalogue manifest. All 23 have no GFF3 / FAA.
  - `old_pipeline_only`: s1–s4 results only.
- **The U–core–D catalogue exists only for Comamonadaceae.** The two other families planned next are
  Microbacteriaceae and Nitrobacteraceae.
- **1,682 vs 1,683:** `manifest.tsv` has 1,683 lines including its header, so the catalogue has **1,682**
  assemblies (1,682 distinct BioSamples and accessions). "1,683" in earlier messages was the line count. In the pool
  the 1,682 split into 1,679 `ucd_checks_done` and 3 `ucd_catalogue_no_pair`.
- Every catalogue assembly's accession.version equals the pool's chosen accession.
