#!/bin/bash
#SBATCH --job-name=pb04_compare
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/kh36969/longread_callout/logs/%x_%A_%a.out
# Within-species comparison groups (plan 3.2): 1 = terricola, 2 = fungorum.
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
HERE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments/paraburkholderia_round1
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
RUNS=$PROJ/runs_paraburkholderia_round1
case $SLURM_ARRAY_TASK_ID in
    1) GROUP=terricola; STRAINS="SMT4a LjRoot223" ;;
    2) GROUP=fungorum; STRAINS="R-93128 R-91497" ;;
esac
args=()
for s in $STRAINS; do
    fa=$(python -c "import json;print(json.load(open('$PROJ/data/assemblies/$s/manifest.json'))['samples'][0]['fasta'])")
    args+=(--strain $s $RUNS/$s/callout $RUNS/$s/evidence $fa)
done
echo "[$(date -Is)] group=$GROUP strains=$STRAINS host=$(hostname)"
nucmer --version
blastn -version | head -1
/usr/bin/time -v python $HERE/compare_strains.py "${args[@]}" --out $RUNS/compare_$GROUP 2> $RUNS/compare_${GROUP}_time.txt
echo "[$(date -Is)] done"
