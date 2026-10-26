#!/bin/bash
#SBATCH --job-name=fna_s5_rc
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%A_%a.out
# v0.3 read check (s3_readcheck.py) on one S4 candidate list of the full S3 set, in three stages:
#   prep       alleles for the list's loci (all calls), then download the linked runs (--jobs concurrent, MD5 checked);
#              writes samples.txt (one sample per line, those with alleles)
#   map        array task i maps sample i of samples.txt to its alleles and to its whole assembly (skips a sample
#              whose genome_<sample>.json exists, so failed tasks can be resubmitted)
#   summarize  verdicts, evidence levels -> repo
# Usage: sbatch s5_readcheck.sbatch <pool> <list: main|partner_mixed|short_other_segment|replacement> <s4 dir> <stage>
#        (for map: sbatch --array=0-<n-1> ...)
set -euo pipefail
POOL=${1:?pool}; LIST=${2:?list}; S4=${3:?s4 dir name}; STAGE=${4:?stage}
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
S3=$PROJ/fna_pool/$POOL/s3_full
W=$PROJ/fna_pool/$POOL/s5_readcheck/$LIST
REPO=$HERE/$POOL/s5_readcheck_$LIST
READS=$PROJ/data/reads
echo "[$(date -Is)] host=$(hostname) stage=$STAGE task=${SLURM_ARRAY_TASK_ID:-} commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
col() { python -c "import csv,sys; print('\n'.join(r['$1'] for r in csv.DictReader(open('$S3/pairs.tsv'), delimiter='\t')))"; }
FASTA=$( (paste -d= <(col sample_a) <(col fna_a); paste -d= <(col sample_b) <(col fna_b)) | sort -u | tr '\n' ' ')
case $STAGE in
prep)
    [ ! -e $W/tests.tsv ] || { echo "Refusing to reuse existing $W" >&2; exit 2; }
    mkdir -p $W $REPO
    cp $HERE/$POOL/$S4/list_$LIST.tsv $W/loci_ids.tsv
    python $HERE/s3_readcheck.py alleles --loci $S3/merged/loci.tsv --calls $S3/merged/calls_merged.tsv --loci-ids $W/loci_ids.tsv \
        --fasta $FASTA --out $W
    ls $W/alleles_*.fa | sed 's/.*alleles_//; s/\.fa$//' > $W/samples.txt
    echo "[$(date -Is)] fetch, $(wc -l < $W/samples.txt) samples"
    python $HERE/s3_readcheck.py fetch --pool-dir $HERE/$POOL --work $W --reads $READS --jobs 6
    ;;
map)
    S=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $W/samples.txt)
    if [ -e $W/paf/genome_$S.json ]; then echo "$S done before"; exit 0; fi
    python $HERE/s3_readcheck.py map --pool-dir $HERE/$POOL --work $W --reads $READS --sample $S --threads $SLURM_CPUS_PER_TASK
    python $HERE/s3_readcheck.py genome --pool-dir $HERE/$POOL --work $W --reads $READS --sample $S --fasta $FASTA --threads $SLURM_CPUS_PER_TASK
    ;;
summarize)
    N=$(wc -l < $W/samples.txt); D=$(ls $W/paf/genome_*.json | wc -l)
    [ "$N" = "$D" ] || { echo "only $D of $N samples mapped" >&2; exit 2; }
    python $HERE/s3_readcheck.py summarize --work $W --out $REPO
    cp $W/paf/*.json $REPO/
    for R in $(ls $W/paf/SAM*.paf 2>/dev/null | sed 's/.*_\([A-Z]*RR[0-9]*\)\.paf/\1/' | sort -u); do cp $READS/$R/verify.json $REPO/verify_$R.json; done
    ;;
esac
echo "[$(date -Is)] done"
