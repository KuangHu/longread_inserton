#!/bin/bash
#SBATCH --job-name=fna_s3_pilot
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%j.out
# S3 small-scale pilot on s3_pilot/pilot_pairs.tsv: call each pair, merge loci, export U/core/D, plot, summarise.
# Usage: sbatch s3_pilot.sbatch <pool name>
set -euo pipefail
POOL=${1:?pool name}
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
REPO=$HERE/$POOL/s3_pilot
W=$PROJ/fna_pool/$POOL/s3_pilot
echo "[$(date -Is)] host=$(hostname) commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
mkdir -p $W
col() { python -c "import csv,sys; print('\n'.join(r['$1'] for r in csv.DictReader(open('$REPO/pilot_pairs.tsv'), delimiter='\t')))"; }
paste <(col pair_id) <(col sample_a) <(col fna_a) <(col sample_b) <(col fna_b) | while IFS=$'\t' read -r PID SA FA SB FB; do
    echo "[$(date -Is)] call $PID $SA $SB"
    /usr/bin/time -v -o $W/time_$PID.txt python $HERE/s3_pair_calls.py call --pair-id $PID --name-a $SA --fasta-a $FA \
        --name-b $SB --fasta-b $FB --out $W/$PID --threads $SLURM_CPUS_PER_TASK
done
FASTA=$(paste -d= <(col sample_a) <(col fna_a); paste -d= <(col sample_b) <(col fna_b))
GFF=$(paste -d= <(col sample_a) <(col gff_a); paste -d= <(col sample_b) <(col gff_b))
FASTA=$(echo "$FASTA" | sort -u | tr '\n' ' '); GFF=$(echo "$GFF" | sort -u | tr '\n' ' ')
echo "[$(date -Is)] merge"
python $HERE/s3_pair_calls.py merge --calls $W/P*/calls.tsv --fasta $FASTA --out $W/merged
cp $W/merged/loci.tsv $W/merged/calls_merged.tsv $REPO/
IDS=$(python $HERE/s3_pilot.py plotlist --loci $W/merged/loci.tsv --calls $W/merged/calls_merged.tsv)
echo "[$(date -Is)] plot $(echo $IDS | wc -w) loci"
python $HERE/s3_pair_calls.py plot --loci $W/merged/loci.tsv --calls $W/merged/calls_merged.tsv --fasta $FASTA --gff $GFF \
    --loci-ids $IDS --out $REPO/plots
cat $W/P*/one_flank.tsv | awk 'NR==1 || !/^pair_id/' > $REPO/one_flank.tsv
python $HERE/s3_pilot.py summary --pilot $REPO/pilot_pairs.tsv --loci $W/merged/loci.tsv --calls $W/merged/calls_merged.tsv \
    --work $W --out $REPO
grep -h "Elapsed\|Maximum res" $W/time_P*.txt | paste - - > $REPO/pair_resources.txt
echo "[$(date -Is)] done"
