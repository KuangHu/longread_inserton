#!/bin/bash
#SBATCH --job-name=fna_s3_readcheck
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%j.out
# S3 pilot long-read allele check on the anchored loci: build alleles, fetch the linked runs, map each sample, summarise.
# Usage: sbatch s3_readcheck.sbatch <pool name>
set -euo pipefail
POOL=${1:?pool name}
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
REPO=$HERE/$POOL/s3_pilot
W=$PROJ/fna_pool/$POOL/s3_pilot/readcheck
READS=$PROJ/data/reads
echo "[$(date -Is)] host=$(hostname) commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
[ ! -e $W/tests.tsv ] || { echo "Refusing to reuse existing $W" >&2; exit 2; }
mkdir -p $W $REPO/readcheck
col() { python -c "import csv,sys; print('\n'.join(r['$1'] for r in csv.DictReader(open('$REPO/pilot_pairs.tsv'), delimiter='\t')))"; }
FASTA=$( (paste -d= <(col sample_a) <(col fna_a); paste -d= <(col sample_b) <(col fna_b)) | sort -u | tr '\n' ' ')
python $HERE/s3_readcheck.py alleles --loci $REPO/loci.tsv --calls $REPO/calls_merged.tsv --fasta $FASTA --out $W
echo "[$(date -Is)] fetch"
python $HERE/s3_readcheck.py fetch --pool-dir $HERE/$POOL --work $W --reads $READS
for FA in $W/alleles_*.fa; do
    S=$(basename $FA .fa); S=${S#alleles_}
    echo "[$(date -Is)] map $S"
    python $HERE/s3_readcheck.py map --pool-dir $HERE/$POOL --work $W --reads $READS --sample $S --threads $SLURM_CPUS_PER_TASK
done
python $HERE/s3_readcheck.py summarize --work $W --out $REPO/readcheck
cp $W/paf/*.json $REPO/readcheck/
for R in $(ls $W/paf/*.paf | sed 's/.*_\([A-Z]*RR[0-9]*\)\.paf/\1/' | sort -u); do cp $READS/$R/verify.json $REPO/readcheck/verify_$R.json; done
echo "[$(date -Is)] done"
