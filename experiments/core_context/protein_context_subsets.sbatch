#!/bin/bash
#SBATCH --job-name=core_ctx_prot
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# Protein entry on the frozen subsets (subsets_20261027): all GFF3 CDS / FAA proteins -> MMseqs families with member
# check -> single- and two-gene cores with oriented neighbour-gene contexts. Work files (proteins.faa, mmseqs) on scratch.
# Usage: sbatch --array=0-2 protein_context_subsets.sbatch
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
rm -rf $W
for step in genes cluster call; do
    extra=(); [ $step = genes ] && extra=(--manifest $M); [ $step != genes ] && extra=(--threads $SLURM_CPUS_PER_TASK)
    [ $step = call ] && extra+=(--out $SEL/$S/protein_entry)
    /usr/bin/time -v -o $W.$step.time python $HERE/protein_context.py $step --work $W "${extra[@]}" > /dev/null
    grep -E "Elapsed|Maximum resident" $W.$step.time
done
echo "[$(date -Is)] done"
