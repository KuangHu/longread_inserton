#!/bin/bash
#SBATCH --job-name=ucd_readcheck
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# Long-read check of U-core-D positions (ucd_readcheck.py), one stage per submission:
#   genome     --array=0-<n-1>  one assembly of <round>/samples.txt per task: whole-run reads -> whole assembly, once
#   check      --array=0-<n-1>  one assembly per task: allele check of the round's positions
#   summarize  one task (FIGURES=<pair ids>)
# Usage: sbatch [--array=..] ucd_readcheck.sbatch <round_dir> <stage> [repo_out]
set -euo pipefail
RD=${1:?round_dir}; STAGE=${2:?stage}; OUT=${3:-}
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/core_context
T=${SLURM_ARRAY_TASK_ID:-0}
S=$(sed -n "$((T + 1))p" $RD/samples.txt)
mkdir -p $RD/time
echo "[$(date -Is)] host=$(hostname) stage=$STAGE task=$T sample=$S commit=$(git -C $CODE rev-parse HEAD) dirty=$(git -C $CODE status --porcelain experiments/core_context | wc -l)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /global/scratch/users/kh36969/longread_callout/env/msi
cd $HERE
TIME="/usr/bin/time -v -o $RD/time/${STAGE}_$T.time"
case $STAGE in
genome)    $TIME python ucd_readcheck.py genome --sample $S --threads $SLURM_CPUS_PER_TASK ;;
check)     $TIME python ucd_readcheck.py check --round-dir $RD --sample $S --threads $SLURM_CPUS_PER_TASK ;;
summarize) $TIME python ucd_readcheck.py summarize --round-dir $RD --out ${OUT:?repo_out} --figure-pairs "${FIGURES:-}" ;;
esac
echo "[$(date -Is)] done"
