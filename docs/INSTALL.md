# Runtime dependencies

Use Linux or WSL, Python 3.10+, a C/C++ compiler, and the native executables below.
No native binaries or third-party source archives are bundled in this package.

| Executable | Pilot version/source |
|---|---|
| `pals` | PALS 1.0, [official source](https://www.drive5.com/pals/pals_source.tar.gz) |
| `piler2` | PILER 1.0, [official site](https://www.drive5.com/piler/) |
| `muscle` | MUSCLE 3.8.1551 (the v3 `-in/-out/-diags1` interface) |
| `RepeatScout`, `build_lmer_table` | Dfam-maintained 1.0.7, commit `a1dcd4022e70fedaa2e28d8d08d4eac898e79738` |
| `blastn`, `blastx`, `makeblastdb` | NCBI BLAST+ 2.17.0 |
| Python `Bio` | Biopython ≥1.86,<2 |

Do not silently substitute MUSCLE 5 or an old RepeatScout lacking `-ranges`.
PALS/PILER are old programs, so successful installation alone is not proof that
they can process an arbitrary input size or repeat density. Failures are logged
and stop a run; the workflow does not replace a failed PALS search with BLAST.

## Maintained RepeatScout

```bash
git clone https://github.com/Dfam-consortium/RepeatScout.git env/RepeatScout
git -C env/RepeatScout checkout a1dcd4022e70fedaa2e28d8d08d4eac898e79738
make -C env/RepeatScout
export PATH="$PWD/env/RepeatScout:$PATH"
```

PALS and PILER sources are linked from their official sites. In the pilot, PALS
was built without changing algorithm source, replacing obsolete CPU/static-link
flags with:

```bash
make -j2 CFLAGS='-O3 -funroll-loops -DNDEBUG=1' LDLIBS=-lm
```

Use the source-specific makefile for PILER; the executable name required here is
`piler2`. Consult the [PILER manual](https://www.drive5.com/piler/piler_userguide.html)
for its build and native DF workflow. Put both executables on PATH.

For BLAST+, use the [official 2.17.0 archive](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/)
for your platform. MUSCLE 3.8.1551 is available through Bioconda; pin the major
version rather than installing the latest MUSCLE interface.

Executable paths can be supplied explicitly:

```bash
multisite-insertion run --manifest samples.json --out runs/example \
  --pals /path/pals --piler2 /path/piler2 --muscle /path/muscle \
  --RepeatScout /path/RepeatScout --build_lmer_table /path/build_lmer_table \
  --blastn /path/blastn --blastx /path/blastx --makeblastdb /path/makeblastdb
```

`doctor` prints path, SHA256 and the tool's version/help banner. Presence is not
a compatibility guarantee; exercise `examples/make_demo.py` after setup. Every
actual command and used executable hash is recorded in each run's provenance.
When both native outputs are imported, only BLAST+ is needed for a full annotated
analysis; no native detector binary is invoked.

## Build a wheel

```bash
python -m pip install build
python -m build
python -m pip install dist/multisite_insertion-0.1.0-py3-none-any.whl
```

Data, environments, and native binaries are excluded from package builds.
