#!/bin/bash
#SBATCH --job-name=core_ctx_scale
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out
# Assembly-only scale-up of the protein entry (scale_entry.py), one stage per submission:
#   genes        --array=0-<n-1>   per-assembly genes, n tasks over the batch (resumable)
#   pool_proteins one task: distinct proteins of every manifest assembly (read-only)
#   families     one task, more CPUs / memory (--cpus-per-task=40 --mem=180G): the pool mapping, clustered once
#   occurrences  one task
#   compare      --array=0-63      one shard per task over every batch in BATCHES (resumable)
#   summarize    one task
#   windows      one task per batch: S6 near / host windows of target-scope occurrences (scale_review.py)
#   s6           --array=0-63      S6 review per shard over BATCHES (needs compare + windows)
#   boundary     --array=0-63      boundary review of the structural and cross-assembly pairs per shard (needs s6);
#                                  REUSE=<tags> reuses pair rows of earlier boundary runs (e.g. REUSE=b1_v1)
#   s6_summarize one task
#   candidates   one task: evidence tiers, core list, validation batch from local reads (scale_candidates.py)
# Usage: sbatch [--array=..] [--cpus-per-task=.. --mem=..] scale_batch.sbatch <work> <manifest> <batch> <batches> <stage> [repo_out]
set -euo pipefail
WORK=${1:?work}; MAN=${2:?manifest}; B=${3:?batch}; BS=${4:?batches}; STAGE=${5:?stage}; OUT=${6:-}
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/core_context
T=${SLURM_ARRAY_TASK_ID:-0}
mkdir -p $WORK/time
echo "[$(date -Is)] host=$(hostname) stage=$STAGE task=$T commit=$(git -C $CODE rev-parse HEAD) dirty=$(git -C $CODE status --porcelain experiments/core_context | wc -l)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /global/scratch/users/kh36969/longread_callout/env/msi
cd $HERE
TIME="/usr/bin/time -v -o $WORK/time/${STAGE}_b${B}_${BS//,/-}_$T.time"
case $STAGE in
genes)       $TIME python scale_entry.py genes --manifest $MAN --work $WORK --batch $B --task $T --ntasks ${SLURM_ARRAY_TASK_COUNT:-1} ;;
pool_proteins) $TIME python scale_entry.py pool_proteins --manifest $MAN --work $WORK ;;
families)    $TIME python scale_entry.py families --pool --manifest $MAN --work $WORK --batch $B --threads $SLURM_CPUS_PER_TASK ;;
occurrences) $TIME python scale_entry.py occurrences --manifest $MAN --work $WORK --batch $B ;;
compare)     $TIME python scale_entry.py compare --work $WORK --batches $BS --shard $T ${MAX_SIG_PAIRS:+--max-sig-pairs $MAX_SIG_PAIRS} ;;
summarize)   $TIME python scale_entry.py summarize --work $WORK --batches $BS --repo-out $OUT --regression $HERE/subsets_20261027/ucd_candidates_v0.1/candidates.tsv ;;
windows)     $TIME python scale_review.py windows --work $WORK --batch $B ;;
s6)          $TIME python scale_review.py s6 --work $WORK --batches $BS --shard $T ;;
boundary)    $TIME python scale_review.py boundary --work $WORK --batches $BS --shard $T --manifest $MAN ${REUSE:+--reuse $REUSE} ;;
s6_summarize) $TIME python scale_review.py summarize --work $WORK --batches $BS --repo-out $OUT --regression $HERE/subsets_20261027/ucd_candidates_v0.1/candidates.tsv ;;
candidates)  $TIME python scale_candidates.py --work $WORK --batches $BS --repo-out $OUT ;;
esac
grep -E "Elapsed|Maximum resident" $WORK/time/${STAGE}_b${B}_${BS//,/-}_$T.time
echo "[$(date -Is)] done"
