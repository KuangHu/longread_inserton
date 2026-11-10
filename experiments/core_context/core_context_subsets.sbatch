#!/bin/bash
#SBATCH --job-name=core_ctx
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# Graph entry on frozen comparable subsets (subsets_20261027): PanGraph --verify build, then repeated-block U/core/D call.
# Array task i = i-th subset in sorted order. Big files (joint FNA, graph JSON, contexts.fna) stay on scratch.
# Usage: sbatch --array=0-2 core_context_subsets.sbatch
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/core_context
SEL=$HERE/subsets_20261027
PG=$PROJ/env/pangraph-1.4.0/pangraph-x86_64-unknown-linux-gnu
M=$(ls $SEL/*.samples.tsv | sort | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
S=$(basename $M .samples.tsv)
W=$PROJ/core_context/subsets_20261027/$S
echo "[$(date -Is)] host=$(hostname) subset=$S commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
mkdir -p $W
[ -e $W/bundle/graph_provenance.json ] || { rm -rf $W/bundle
    /usr/bin/time -v -o $W/build.time python $HERE/core_context.py build --manifest $M --out $W/bundle --pangraph $PG --threads $SLURM_CPUS_PER_TASK; }
rm -rf $W/calls
/usr/bin/time -v -o $W/call.time python $HERE/core_context.py call --bundle $W/bundle --out $W/calls > /dev/null
grep -E "Elapsed|Maximum resident" $W/build.time $W/call.time
echo "[$(date -Is)] done"
