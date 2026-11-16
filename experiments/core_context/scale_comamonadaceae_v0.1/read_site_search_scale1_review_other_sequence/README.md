# Scale 1 supplementary review: `not_in_assembly:other_sequence`

The scale-1 batch has one not-in-assembly candidate site of class `not_in_assembly:other_sequence`:
SAMN42749333 PF_7dd788787bef, NZ_CP162616.1:456,191–458,920, 3 molecules.
- The reads carry the core between two flanks.
- The assembly holds about 2.7 kb of different sequence there.

`candidate_evidence.py` reviewed `core_absent` only. This folder runs the same evidence and full-insertion steps,
with thresholds unchanged, using `--classes not_in_assembly:other_sequence`. The alleles are built the same way:
- **P:** left flank + the reads' insert + right flank.
- **A:** the assembly as it is, which here includes the 2.7 kb segment.

Inputs:
- `positions.tsv.gz` and `samples.tsv` are the scale-1 rows of this sample.
- The work directory links the scale-1 per-sample evidence (`samples/SAMN42749333`) and the panel, read-only.

The scale-1 batch results are not changed.
