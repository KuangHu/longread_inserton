#!/bin/bash
#SBATCH --job-name=core_ctx_s6rev
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# S6 v0.1 review layer (s6_review.py): boundary denominators, same-assembly context classes, unit relations, cases.
# Usage: sbatch --array=0-2 --exclude=<list> s6_review_subsets.sbatch
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/core_context
SEL=$HERE/subsets_20261027
M=$(ls $SEL/*.samples.tsv | sort | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
S=$(basename $M .samples.tsv)
W=$PROJ/core_context/subsets_20261027/$S/protein
echo "[$(date -Is)] host=$(hostname) subset=$S commit=$(git -C $CODE rev-parse HEAD) dirty=$(git -C $CODE status --porcelain experiments/core_context | wc -l)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
cd $HERE
/usr/bin/time -v -o $W.s6review.time python s6_review.py --subset-dir $SEL/$S --work $W > /dev/null
grep -E "Elapsed|Maximum resident" $W.s6review.time
echo "[$(date -Is)] done"
