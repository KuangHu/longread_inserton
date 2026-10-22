# Changelog

## 0.1.0 — native multisite repeat route

- New `multisite-insertion` distribution and CLI; independent `multisite_insertion` Python module.
- Direct default PALS/PILER-DF and maintained RepeatScout runs, plus immutable native-output imports.
- Post-hoc deposited annotation, tnp-supported/all subsets, candidate boundary inference and flank exports.
- Separate native, inferred and reference-annotation coordinates; retain unresolved, censored and invalid records.
- Command/input/executable hashes, sequence validation, regression tests and a deterministic synthetic example.
- Replaces the former `insertion-site-recovery` package in this repository. No old reads/SV,
  empty/filled, Arm A or Arm B implementation is part of the new distribution.
- Old code and results remain reachable at commit `5ad842b0c8261485c27fa685de9fc538ff49cd23`.
