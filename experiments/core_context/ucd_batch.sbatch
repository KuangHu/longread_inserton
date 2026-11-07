#!/bin/bash
#SBATCH --job-name=ucd_batch
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# usage: sbatch [--array=0-47] ucd_batch.sbatch <freeze|chunk|merge> [extra args]
# OUT/WORK default to the b1-2-3 match_both_lt2500 run; a resubmitted array skips finished chunks.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /global/scratch/users/kh36969/longread_callout/env/msi
cd /global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments/core_context
OUT=${OUT:-scale_comamonadaceae_v0.1/ucd_batch_v0.1}
WORK=${WORK:-/global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.1}
export TMPDIR=${SLURM_TMPDIR:-/tmp}
CMD=$1; shift
/usr/bin/time -v python ucd_batch.py $CMD --out $OUT --work $WORK "$@"
