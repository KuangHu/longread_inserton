#!/bin/bash
#SBATCH --job-name=dna_ucd
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --exclude=n0004.lr6,n0065.lr6,n0072.lr6,n0018.lr6,n0054.lr6
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# DNA-first discovery entry (dna_discovery.py), one stage per submission. Reads verified FNA + manifest only.
#   seeds        --array=0-<chunks-1>  --cpus-per-task=32 --mem=90G   all-vs-all per seed chunk
#   library      one task                                              new library shard
#   occurrences  --array=<dna_discovery.py occurrences --list>         (library shard x genome part)
#   pairs        --array=0-<64*SUBS-1>                                 (v0.1.2 legacy) extension per library-sequence shard
#   enumerate    --array=0-63  --gen G [--max-pairs N]                 v0.2 source pairs -> alignment-key shards
#   align        --array=0-127 --gen G --cpus-per-task=16              v0.2 each distinct key once, then extension
# Usage: sbatch [..] dna_batch.sbatch <work> <manifest> <stage> [extra args]
set -euo pipefail
WORK=${1:?work}; MAN=${2:?manifest}; STAGE=${3:?stage}; shift 3
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/core_context
T=$(( ${SLURM_ARRAY_TASK_ID:-0} + ${OFFSET:-0} ))  # OFFSET: arrays above MaxArraySize (1001) are submitted in slices
mkdir -p $WORK/time
echo "[$(date -Is)] host=$(hostname) stage=$STAGE task=$T commit=$(git -C $CODE rev-parse HEAD) dirty=$(git -C $CODE status --porcelain experiments/core_context | wc -l)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /global/scratch/users/kh36969/longread_callout/env/msi
export TMPDIR=${TMPDIR:-/tmp}
cd $HERE
# PACK=n: array task T runs the units T*n .. T*n+n-1 side by side on the node, each with cpus/n threads (lr6 nodes are
# exclusive: a job is charged the whole node, so small units share one job). Finished units skip themselves.
PACK=${PACK:-1}
TH=$(( SLURM_CPUS_PER_TASK / PACK )); [ $TH -ge 1 ] || TH=1
pids=()
for i in $(seq 0 $((PACK - 1))); do
  U=$(( T * PACK + i ))
  /usr/bin/time -v -o $WORK/time/${STAGE}_${SLURM_JOB_ID}_$U.time \
    python dna_discovery.py $STAGE --work $WORK --manifest $MAN --index $U --threads $TH "$@" &
  pids+=($!)
done
rc=0; for p in "${pids[@]}"; do wait $p || rc=1; done
grep -hE "Elapsed|Maximum resident" $WORK/time/${STAGE}_${SLURM_JOB_ID}_*.time | sort | tail -2
[ $rc -eq 0 ] || { echo "[$(date -Is)] a unit failed"; exit 1; }
echo "[$(date -Is)] done"
