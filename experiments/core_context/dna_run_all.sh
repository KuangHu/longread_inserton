#!/bin/bash
# All pools through the DNA entry, one batch after another in one work dir (each batch searched against the
# earlier ones). Order: Microbacteriaceae (4), Nitrobacteraceae (5), Comamonadaceae (1-3), then 6, 7, 8.
# Stops at the first batch that does not finish; rerunning resumes (every stage skips finished units).
# usage: dna_run_all.sh <work> <manifest> [max_pairs_per_bucket] [max_gens]
WORK=${1:?work}; MAN=${2:?manifest}; MAXP=${3:-4000000}; MAXG=${4:-3}
HERE=$(cd $(dirname $0) && pwd)
for spec in b4:4 b5:5 b1:1,2,3 b6:6 b7:7 b8:8; do
  tag=${spec%%:*}; batches=${spec#*:}
  if [ -f $WORK/report_$tag.json ]; then echo "[$(date -Is)] $tag already done"; continue; fi
  echo "[$(date -Is)] start $tag (batches $batches)"
  bash $HERE/dna_run.sh $WORK $MAN $tag $batches $MAXP $MAXG || { echo "[$(date -Is)] STOP at $tag"; exit 1; }
done
echo "[$(date -Is)] all batches done"
