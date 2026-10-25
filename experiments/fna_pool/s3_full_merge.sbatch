#!/bin/bash
#SBATCH --job-name=fna_s3_merge
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%j.out
# Full S3 candidate discovery, merge: all pairs' raw calls -> loci (difference_class with --fasta), per-pair summary,
# short-locus review table. Large tables stay on scratch; the repo gets summaries, the short anchored review table
# and sha256 of the full tables. Usage: sbatch s3_full_merge.sbatch <pool name>
set -euo pipefail
POOL=${1:?pool name}
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
W=$PROJ/fna_pool/$POOL/s3_full
REPO=$HERE/$POOL/s3_full
echo "[$(date -Is)] host=$(hostname) commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
mkdir -p $REPO
N=$(($(wc -l < $W/pairs.tsv) - 1)); DONE=$(ls $W/pairs/*/call_run.json | wc -l)
echo "pairs $N, called $DONE"
[ "$N" = "$DONE" ] || { echo "Not all pairs called" >&2; exit 2; }
col() { python -c "import csv,sys; print('\n'.join(r['$1'] for r in csv.DictReader(open('$W/pairs.tsv'), delimiter='\t')))"; }
FASTA=$( (paste -d= <(col sample_a) <(col fna_a); paste -d= <(col sample_b) <(col fna_b)) | sort -u | tr '\n' ' ')
echo "[$(date -Is)] merge"
python $HERE/s3_pair_calls.py merge --calls $W/pairs/P*/calls.tsv --fasta $FASTA --out $W/merged
cat $W/pairs/P*/one_flank.tsv | awk 'NR==1 || !/^pair_id/' > $W/merged/one_flank.tsv
mkdir -p $W/work && for P in $W/pairs/P*; do ln -sfn $P $W/work/$(basename $P); done
python $HERE/s3_pilot.py summary --pilot $W/pairs.tsv --loci $W/merged/loci.tsv --calls $W/merged/calls_merged.tsv \
    --work $W/work --out $W/merged
python $HERE/s3_pilot.py crosstab --pilot $W/pairs.tsv --loci $W/merged/loci.tsv --out $W/merged
cp $W/merged/pilot_summary.tsv $REPO/pair_summary.tsv
cp $W/merged/pilot_summary.json $REPO/summary.json
cp $W/merged/short_loci_review.tsv $W/pairs.tsv $REPO/
cat $W/pairs/P*/call_run.json | python -c "
import json, sys
from collections import Counter
runs = [json.loads(x) for x in sys.stdin.read().replace('}\n{', '}\x00{').split('\x00')]
print(json.dumps(dict(pairs=len(runs), callers=dict(Counter(r['caller'] for r in runs)), minimap2=dict(Counter(r['minimap2'] for r in runs))), indent=1))
" > $REPO/call_runs.json
(cd $W/merged && sha256sum loci.tsv calls_merged.tsv one_flank.tsv ucd.fasta) > $REPO/full_tables.sha256
grep -h "Elapsed\|Maximum res" $W/time/P*.txt | paste - - > $REPO/pair_resources.txt
echo "[$(date -Is)] done"
