#!/bin/bash
#SBATCH --job-name=d_fetch
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/global/home/users/kh36969/logs/core_context/%x_%j.out
# Fetch job of one family for the DNA entry: catalog -> pool -> download + verification -> manifest -> queue line.
# The only job that downloads. It holds ONE CPU on the shared partition lr7 (charged one CPU; lr6 is node-exclusive
# and would be charged a whole node), at most 24 h, i.e. at most 24 CPU-hours per submission. download_fna.py
# refuses to run in SLURM unless FETCH_JOB=1 and the job holds exactly one CPU. Analysis jobs never download: they
# read the manifest written here, whose FASTA hashes are re-checked by the run.
# Resumable: resubmitting continues the download (verified accessions are not fetched again).
#
# usage: sbatch dna_fetch.sbatch <family_name> <ncbi_taxon_id> <fetch_root> [queue.tsv]
#   <fetch_root>/<family>/{raw,pool,data,manifest.tsv}; with queue.tsv, a line for dna_queue.sh is appended once
#   LIMIT=n (env): download only the first n accessions (pilot); the queue line is then not written
#   CAP=n (env): keep at most n BioSamples of the pool, chosen by dna_cap.py (metadata-only rule dna_cap_v1)
set -euo pipefail
FAM=${1:?family}; TAXON=${2:?taxon}; ROOT=${3:?fetch root}; QUEUE=${4:-}
[ "${SLURM_CPUS_ON_NODE:-}" = 1 ] || { echo "fetch job must hold exactly one CPU (got ${SLURM_CPUS_ON_NODE:-none})"; exit 2; }
export FETCH_JOB=1
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments
D=$ROOT/$FAM; mkdir -p $D/pool
say() { echo "[$(date -Is)] $*"; }
say "fetch $FAM taxon=$TAXON host=$(hostname) cpus=$SLURM_CPUS_ON_NODE commit=$(git -C $CODE rev-parse --short HEAD)"
cd $CODE/fna_pool
if [ ! -f $D/pool/pool.tsv ]; then
  /usr/bin/python3 build_pool.py fetch --taxon $TAXON --work $D
  /usr/bin/python3 build_pool.py build --taxon $TAXON --work $D --out $D/pool --selection assembly_first_v1
fi
# the DNA entry needs one assembly per BioSample: only the primary download rows (no GCA identity copies)
POOL=$D/pool/pool.tsv
if [ -n "${CAP:-}" ]; then
  POOL=$D/pool/pool_cap$CAP.tsv
  [ -f $POOL ] || /usr/bin/python3 $CODE/core_context/dna_cap.py --pool $D/pool/pool.tsv --cap $CAP --out $POOL
  [ -f $D/pool/download_primary.tsv ] ||
    awk -F'\t' 'NR==FNR{if(FNR>1)k[$3]=1; next} FNR==1 || ($3=="primary" && ($1 in k))' $POOL $D/pool/download_manifest.tsv > $D/pool/download_primary.tsv
fi
[ -f $D/pool/download_primary.tsv ] ||
  awk -F'\t' 'NR==1 || $3=="primary"' $D/pool/download_manifest.tsv > $D/pool/download_primary.tsv
say "pool: $(($(wc -l < $D/pool/pool.tsv) - 1)) BioSamples, $(($(wc -l < $D/pool/download_primary.tsv) - 1)) accessions to hold"
LIM=""; [ -z "${LIMIT:-}" ] || LIM="--limit $LIMIT"
/usr/bin/python3 download_fna.py run --manifest $D/pool/download_primary.tsv --data $D/data --workers 4 $LIM 2>&1 | tee -a $D/download.log | grep -E "^\[" || true
cd $CODE/core_context
/usr/bin/python3 dna_manifest.py --pool $POOL --data $D/data --family $FAM --out $D/manifest.tsv
N=$(($(wc -l < $D/manifest.tsv) - 1)); P=$(($(wc -l < $POOL) - 1))
say "manifest: $N of $P pool BioSamples verified ($(du -sh $D/data | cut -f1))"
if [ -n "$QUEUE" ] && [ -z "${LIMIT:-}" ] && [ $N -gt 0 ]; then
  grep -q "^$FAM	" $QUEUE 2>/dev/null || { printf '%s\t1\t%s\n' $FAM $D/manifest.tsv >> $QUEUE; say "queued $FAM in $QUEUE"; }
fi
say "done"
