# Batch-1 S6 review, strict rule v0 (version record)

Produced by commit `4175584` (`scale_review.py s6 / boundary / summarize`, batch 1). Kept unchanged as a record.

In this version a same-assembly pair was a U–core–D pair candidate only if both copies were S6 supported positions
(`site_group` or `single_observation_site`). It also had to have host context different on both sides and no host gene
shared physically. That position requirement was later removed as a hard gate. S6 grouping is now a separate evidence
field (see `../README.md`), because cross-assembly grouping status changes as assemblies are added, while the pair's own
structure does not.
