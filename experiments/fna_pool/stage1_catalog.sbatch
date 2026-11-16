#!/bin/bash
#SBATCH --job-name=fna_stage1_catalog
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%j.out
# Section G: cross-batch protein dedup -> global MMseqs2 families with member checks -> catalog and stage-2 manifests.
# Usage: sbatch stage1_catalog.sbatch
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
W=$PROJ/fna_pool/scaleup_stage1_20261026/stage1_catalog
COMMIT=$(git -C $CODE rev-parse HEAD)
echo "[$(date -Is)] host=$(hostname) commit=$COMMIT"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
[ ! -e $W/mmseqs/cluster_run.json ] || { echo "Refusing to reuse existing $W" >&2; exit 2; }
python $HERE/stage1_catalog.py proteins --work $W
echo "[$(date -Is)] cluster"
python $HERE/stage1_catalog.py cluster --work $W --threads $SLURM_CPUS_PER_TASK
echo "[$(date -Is)] catalog"
python $HERE/stage1_catalog.py catalog --work $W --commit $COMMIT
echo "[$(date -Is)] done"
