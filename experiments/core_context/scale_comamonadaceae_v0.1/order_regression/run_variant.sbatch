#!/bin/bash
#SBATCH --job-name=ucd_order_regr
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# One variant of make_variants.py per array task: check stage of SAMN31258398 with that positions.tsv.
# Usage: sbatch --array=0-<n-1> run_variant.sbatch <regression_dir>
set -euo pipefail
D=${1:?regression_dir}
V=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $D/variants.txt)
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
echo "[$(date -Is)] host=$(hostname) variant=$V commit=$(git -C $CODE rev-parse HEAD) dirty=$(git -C $CODE status --porcelain experiments/core_context | wc -l)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /global/scratch/users/kh36969/longread_callout/env/msi
cd $CODE/experiments/core_context
python ucd_readcheck.py check --round-dir $D/$V --sample SAMN31258398 --threads $SLURM_CPUS_PER_TASK
echo "[$(date -Is)] done"
