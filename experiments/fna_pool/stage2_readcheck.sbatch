#!/bin/bash
#SBATCH --job-name=fna_stage2_rc
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%A_%a.out
# Stage 2 batch 1: known-site read check (s3_readcheck v0.3) of stage-1 main loci with local verified carrier reads.
# No downloads: there is no fetch stage, and only runs re-checked on disk (runs_allowed.tsv) are mapped.
#   prep       select targets from the stage-1 catalog, build alleles for all calls of those loci, re-check local runs
#   map        array task i: sample i of samples.txt -> alleles and whole assembly (skips samples already done)
#   summarize  v0.3 verdicts with per-molecule support records, then the read-evidence layer, comparison and next batch
#   diagnose   <sample>: whole run mapped to the sample's own assembly, facts on why no read spans its sites
#   followup   open items of this batch, read types behind each side, priority-1 download chunks (metadata only)
#   regress    summarize the baseline pilot work dir again with this code; outputs must equal the committed pilot
# Usage: sbatch stage2_readcheck.sbatch <stage> [sample]   (map: sbatch --array=0-<n-1> ...)
set -euo pipefail
STAGE=${1:?stage}
B=comamonadaceae
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
POOLDIR=$HERE/scaleup_stage1_20261026/$B
S3=$PROJ/fna_pool/scaleup_stage1_20261026/$B/s3_full
W=$PROJ/fna_pool/scaleup_stage1_20261026/stage2_readcheck_batch1
REPO=$HERE/scaleup_stage1_20261026/stage2_readcheck_batch1
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
    python $HERE/stage2_readcheck.py select --work $W --batch $B --catalog-commit 23f195b338128f0dd4a39a26b4c75e0f93a1c658
    python $HERE/s3_readcheck.py alleles --loci $S3/merged/loci.tsv --calls $S3/merged/calls_merged.tsv --loci-ids $W/targets.tsv \
        --fasta $FASTA --out $W
    ls $W/alleles_*.fa | sed 's/.*alleles_//; s/\.fa$//' > $W/samples.txt
    python $HERE/stage2_readcheck.py localcheck --work $W --pool-dir $POOLDIR --reads $READS --jobs $SLURM_CPUS_PER_TASK
    ;;
map)
    S=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $W/samples.txt)
    if [ -e $W/paf/genome_$S.json ]; then echo "$S done before"; exit 0; fi
    python $HERE/s3_readcheck.py map --pool-dir $POOLDIR --work $W --reads $READS --sample $S --runs-allowed $W/runs_allowed.tsv \
        --threads $SLURM_CPUS_PER_TASK
    python $HERE/s3_readcheck.py genome --pool-dir $POOLDIR --work $W --reads $READS --sample $S --runs-allowed $W/runs_allowed.tsv \
        --fasta $FASTA --threads $SLURM_CPUS_PER_TASK
    ;;
summarize)
    N=$(wc -l < $W/samples.txt); D=$(ls $W/paf/genome_*.json | wc -l)
    [ "$N" = "$D" ] || { echo "only $D of $N samples mapped" >&2; exit 2; }
    python $HERE/s3_readcheck.py summarize --work $W --out $REPO --records
    python $HERE/stage2_readcheck.py layer --work $W --out $REPO --batch $B
    cp $W/targets.tsv $W/targets.json $W/tests.tsv $W/runs_checked.tsv $W/runs_checked.json $W/runs_allowed.tsv $REPO/
    mkdir -p $REPO/paf_json && cp $W/paf/*.json $REPO/paf_json/
    ;;
diagnose)
    python $HERE/stage2_readcheck.py diagnose --work $W --out $REPO --reads $READS --sample ${2:?sample} --jobs $SLURM_CPUS_PER_TASK
    ;;
followup)
    python $HERE/stage2_readcheck.py followup --work $W --out $REPO
    ;;
regress)
    BW=$PROJ/fna_pool/comamonadaceae_20261024/s5_readcheck/main
    T=$PROJ/fna_pool/scaleup_stage1_20261026/stage2_regress_v03; rm -rf $T; mkdir -p $T
    python $HERE/s3_readcheck.py summarize --work $BW --out $T --records > /dev/null
    for f in readcheck_tests.tsv readcheck_loci.tsv readcheck_summary.json; do
        cmp -s $T/$f $HERE/comamonadaceae_20261024/s5_readcheck_main/$f && echo "$f identical" || echo "$f DIFFERS"
    done
    zcat $T/support_records.tsv.gz | wc -l
    ;;
esac
echo "[$(date -Is)] done"
