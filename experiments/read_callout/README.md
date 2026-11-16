# read_callout experiment

Assembly-only callout plus long-read evidence for four isolates with matched public reads.
Results and interpretation: [docs/READ_CALLOUT_RESULTS.md](../../docs/READ_CALLOUT_RESULTS.md).
This is an experiment, not part of the `multisite-insertion` package.

| File | Role |
|---|---|
| `samples.tsv` | sample, ENA run, GCA/GCF assembly, platform, minimap2 preset, adapter splitting |
| `common.sh` | paths, conda env and per-task sample row (edit `PROJ`/`CODE` for another cluster) |
| `00_prepare.sh` | **transfer node (lrc-xfer) only:** reads via `scripts/reads_ready.py prepare`; GCA+GCF packages via NCBI `datasets` |
| `01_inputs.sbatch` | input check (ready manifest, assembly packages present; `input_not_ready`, exit 3, otherwise); read-unit diagnosis; sequence identity check; package manifest |
| `02_callout.sbatch` | `multisite-insertion run --subset all` on the GCF assembly |
| `03_reads.sbatch` | molecule-named reads → minimap2 (`-c --secondary=yes -N 50 -p 0.8`) → `classify_reads.py` |
| `molecule_reads.py` | rename reads to `<molecule>/<n>`; split PacBio records at internal SMRTbell adapters; adapter diagnostics |
| `classify_reads.py` | per-candidate molecule evidence, non-reference insertions, depth curves |
| `summarize.py` | cross-sample `comparison.tsv` / `comparison.json` |
| `results/` | small per-sample outputs copied from the cluster run (paths shown as `$PROJ`) |

Run order (LBL Lawrencium, account `pc_rubinlab`, partition `lr6`):

```bash
# on lrc-xfer (OTP login, in tmux; rerun to resume), never inside SLURM:
bash 00_prepare.sh
# then on the cluster; no step downloads anything:
J1=$(sbatch --parsable 01_inputs.sbatch)
J2=$(sbatch --parsable --dependency=aftercorr:$J1 02_callout.sbatch)
sbatch --dependency=aftercorr:$J2 03_reads.sbatch
python summarize.py
```

The environment needs the package's native tools (see `docs/INSTALL.md`) plus minimap2 2.30,
samtools and NCBI `datasets`. Raw reads, PAF and intermediate FASTQ stay on scratch.
`classify_reads.py` thresholds are module constants fixed before the first run; change them
only as a new, separately reported version.
