#!/bin/bash
#SBATCH --job-name=read_site_search
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# Read-based site search pilot (read_site_search.py), one stage per submission:
#   panel    one task
#   search   --array=0-<n-1>  one pilot sample (read_site_search_pilot/samples.tsv) per task
#   combine  one task
# Usage: sbatch [--array=..] read_site_search.sbatch <stage>
set -euo pipefail
STAGE=${1:?stage}
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/core_context
T=${SLURM_ARRAY_TASK_ID:-0}
S=$(awk -F'\t' -v t=$((T + 2)) 'NR == t {print $1}' $HERE/scale_comamonadaceae_v0.1/read_site_search_pilot/samples.tsv)
echo "[$(date -Is)] host=$(hostname) stage=$STAGE task=$T sample=$S commit=$(git -C $CODE rev-parse HEAD) dirty=$(git -C $CODE status --porcelain experiments/core_context | wc -l)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /global/scratch/users/kh36969/longread_callout/env/msi
cd $HERE
case $STAGE in
panel)   python read_site_search.py panel ;;
search)  python read_site_search.py search --sample $S --threads $SLURM_CPUS_PER_TASK ;;
combine) python read_site_search.py combine ;;
esac
echo "[$(date -Is)] done"
