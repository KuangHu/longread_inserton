#!/bin/bash
# Transfer node only (lrc-xfer, outside SLURM): reads and assembly of the HiFi pilot, so that pilot.sbatch only
# checks inputs and computes. Reads: scripts/reads_ready.py prepare; assembly: prepare_assembly.py + verify_assembly.py
# (NCBI datasets), kept when manifest.json already exists. Rerun the same command to resume.
set -euo pipefail
[ -z "${SLURM_JOB_ID:-}" ] || { echo "prepare_inputs refused inside SLURM job $SLURM_JOB_ID: run it on lrc-xfer" >&2; exit 2; }
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/pacbio_hifi_pilot
SAMPLE=Sam8H1; RUN=SRR13837247; GCA=GCA_016904925.1; GCF=GCF_016904925.1
BIOSAMPLE=SAMN17825431; TOTAL_BP=4131713
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
python $CODE/scripts/reads_ready.py prepare --runs $RUN
ASM=$PROJ/data/assemblies_pacbio_pilot/$SAMPLE
[ -f $ASM/manifest.json ] || python $CODE/experiments/callout_scale_v1/batch10/prepare_assembly.py --sample $SAMPLE --gca $GCA --gcf $GCF --out $ASM
[ -f $ASM/assembly_verify.json ] || python $HERE/verify_assembly.py --asm-dir $ASM --biosample $BIOSAMPLE --total-bp $TOTAL_BP
