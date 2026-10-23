# Shared settings for the read_callout SLURM steps; sourced by each sbatch script.
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/read_callout
DATASETS=/global/home/users/kh36969/.conda/envs/claude-env/bin/datasets
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
export PATH=$PROJ/env/bin:$PATH

# Row N (1-based, after header) of samples.tsv -> shell variables.
row=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' $HERE/samples.tsv)
IFS=$'\t' read -r SAMPLE RUN GCA GCF PLATFORM PRESET SPLIT <<< "$row"
ASM_DIR=$PROJ/data/assemblies/$SAMPLE
READ_DIR=$PROJ/data/reads
OUT=$PROJ/runs/$SAMPLE
echo "[$(date -Is)] $SLURM_JOB_NAME task=$SLURM_ARRAY_TASK_ID sample=$SAMPLE run=$RUN host=$(hostname)"
