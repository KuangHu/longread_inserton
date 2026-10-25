#!/bin/bash
#SBATCH --job-name=fna_s2
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%j.out
# S2 for one pool: catalog + joint FNA, skani triangle, pair table. Usage: sbatch s2.sbatch <pool name>
set -euo pipefail
POOL=${1:?pool name, e.g. comamonadaceae_20261024}
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
W=$PROJ/fna_pool/$POOL
echo "[$(date -Is)] host=$(hostname) pool=$POOL commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
ARGS="--work $W --pool-dir $HERE/$POOL --out $HERE/$POOL/s2"
for step in catalog ani pairs; do
    echo "[$(date -Is)] $step"
    /usr/bin/time -v -o $W/s2_time_$step.txt python $HERE/build_s2.py $step $ARGS --threads $SLURM_CPUS_PER_TASK
done
echo "[$(date -Is)] done"
