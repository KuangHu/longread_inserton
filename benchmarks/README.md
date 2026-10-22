# Frozen environmental-isolate regression

The package was compared with the completed native-only pilot and its boundary
analysis, using exactly the same original FASTA and deposited annotation records.
Native outputs were imported; BLASTX annotation, context BLASTN, boundary
inference, flank BLASTN and summaries were rerun by the package.

| Sample | Accessions | PILER all / tnp-supported | RS all / tnp-supported |
|---|---|---:|---:|
| AltDE | CP001103.3 | 6 / 4 | 39 / 13 |
| UT26S | AP010803.1–AP010807.1 | 6 / 5 | 44 / 7 |
| RHA1 | CP000431.1–CP000434.1 | 11 / 8 | 155 / 32 |

Total: 261 native families; 69 tnp-supported families with 300 native records
(72 PILER + 228 RS). Records are not deduplicated loci. All 300 are in bounds.
The complete native discovery table contains additional invalid RS records;
selection does not silently repair or drop their original coordinates.

Validated packaging invariants:

- all 300 original, inferred and exported coordinate/status rows match;
- all eight FASTA exports are byte-identical to the preceding analysis;
- family summaries, all flank pair summaries, and cross-tool overlap tables match;
- every exported sequence reconstructs from its original genome and strand;
- original detector files remain unchanged through post-processing.

`pilot_expected.json` freezes the coordinate rows and sequence SHA256 values,
not just the total count. Check a reproduced run with:

```bash
python benchmarks/check_pilot.py --run /path/to/reproduced_pilot
multisite-insertion validate --run /path/to/reproduced_pilot
```

Use the original `Environmental_isolates_PILER_RepeatScout.zip` archive for the
full native-output regression:

```bash
python benchmarks/make_pilot_manifest.py --pilot-root /path/environmental_pilot --out data/pilot.json
multisite-insertion run --manifest data/pilot.json --out runs/pilot
python benchmarks/check_pilot.py --run runs/pilot
```

The manifest builder verifies the frozen input FASTA and annotation hashes.
`source_manifest.json` records public ENA URLs, versions and original download
hashes. `Tnp_supported_insertion_boundaries.zip` contains the
preceding analysis and its evidence. These large data archives and native
binaries are deliberately not part of the source distribution.

The native discovery steps were additionally exercised from the new CLI on the
deterministic fixture from `examples/make_demo.py`. That fixture uses random
DNA with repeated segments and a deliberately artificial protein annotation;
it is a software integration test, not a biological positive control.

## Biological calibration in the original pilot

The independently inferred PILER boundaries for 12 UT26S IS6100 records differ
from deposited complete-element coordinates by at most 1 bp at either end.
Three ISSp1 records differ by 3–4 bp at either end. Export adopts those annotated
boundaries while preserving the independent estimates. Other candidates are not
certified by this small calibration.

PILER's 17 selected families contain 154 within-family member pairs; 146 have
no significant corresponding-flank HSP on either side under the fixed settings.
One pair has physically overlapping genomic flank windows and is flagged.
This description is not a specificity estimate: native PILER-DF already favors
isolated repeat copies, and no independent negative benchmark was evaluated.

The 300 records comprise 31 annotated ranges, 262 provisional homology ranges,
3 window-censored ranges, and 4 unresolved native intervals. Each status remains
visible; the last seven are not complete insertions. A family-level tnp label
does not certify every copy: 67/72 PILER and 175/228 RS exported intervals fully
contain at least one deposited transposase CDS.

`regression_result.json` records the measured package checks and tool identities.
