#!/bin/bash
#SBATCH --job-name=fna_s3_calls
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --array=0-24
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%A_%a.out
# Full S3 candidate discovery, calls: array task i calls pairs with row index % 25 == i of s3_full/pairs.tsv
# (s3_pilot.py allpairs). A pair with call_run.json is skipped, so a failed task can be resubmitted.
# Usage: sbatch s3_full_calls.sbatch <pool name>
set -euo pipefail
POOL=${1:?pool name}
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
W=$PROJ/fna_pool/$POOL/s3_full
echo "[$(date -Is)] host=$(hostname) task=$SLURM_ARRAY_TASK_ID commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
mkdir -p $W/time $W/pairs
python -c "
import csv
for i, r in enumerate(csv.DictReader(open('$W/pairs.tsv'), delimiter='\t')):
    if i % 25 == $SLURM_ARRAY_TASK_ID:
        print('\t'.join(r[k] for k in ('pair_id', 'sample_a', 'fna_a', 'sample_b', 'fna_b')))
" | while IFS=$'\t' read -r PID SA FA SB FB; do
    if [ -e $W/pairs/$PID/call_run.json ]; then echo "$PID done before"; continue; fi
    echo "[$(date -Is)] call $PID $SA $SB"
    /usr/bin/time -v -o $W/time/$PID.txt python $HERE/s3_pair_calls.py call --pair-id $PID --name-a $SA --fasta-a $FA \
        --name-b $SB --fasta-b $FB --out $W/pairs/$PID --threads $SLURM_CPUS_PER_TASK
    rm -f $W/pairs/$PID/*.paf  # whole-genome PAFs are regenerable; calls.tsv keeps the coordinates
done
echo "[$(date -Is)] done"
