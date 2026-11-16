#!/bin/bash
#SBATCH --job-name=fna_s4
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/global/home/users/kh36969/logs/fna_pool/%x_%j.out
# S4 on one S3 result set: review (annotation, translation vs FAA, boundary relations, lists) -> families (MMseqs2)
# -> seeds. Usage: sbatch s4.sbatch <pool name> <s3 set: s3_pilot | s3_full> [output name, default s4_<set>]
set -euo pipefail
POOL=${1:?pool name}
SET=${2:?s3_pilot or s3_full}
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/fna_pool
if [ $SET = s3_pilot ]; then
    S3=$HERE/$POOL/s3_pilot; PAIRS=$S3/pilot_pairs.tsv; LOCI=$S3/loci.tsv; CALLS=$S3/calls_merged.tsv
else
    S3=$PROJ/fna_pool/$POOL/s3_full; PAIRS=$S3/pairs.tsv; LOCI=$S3/merged/loci.tsv; CALLS=$S3/merged/calls_merged.tsv
fi
NAME=${3:-s4_$SET}
W=$PROJ/fna_pool/$POOL/$NAME
REPO=$HERE/$POOL/$NAME
echo "[$(date -Is)] host=$(hostname) commit=$(git -C $CODE rev-parse HEAD)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
[ ! -e $W/s4_loci.tsv.gz ] || { echo "Refusing to reuse existing $W" >&2; exit 2; }
mkdir -p $W $REPO
python $HERE/s4_annotate.py review --loci $LOCI --calls $CALLS --pairs $PAIRS --out $W
echo "[$(date -Is)] families"
python $HERE/s4_annotate.py families --pairs $PAIRS --s4 $W --work $W/mmseqs --mmseqs $PROJ/env/s4/bin/mmseqs \
    --threads $SLURM_CPUS_PER_TASK --out $W
echo "[$(date -Is)] seeds"
python $HERE/s4_annotate.py seeds --pairs $PAIRS --s4 $W --work $W --out $W
cp $W/s4_loci.tsv.gz $W/annotation_files.tsv $W/s4_summary.json $W/list_*.tsv $W/families.tsv $W/families_summary.json \
   $W/candidate_loci_families.tsv $W/protein_members.tsv $W/seed_catalog.tsv $REPO/
(cd $W && sha256sum seed_catalog.fa mmseqs/proteins.faa) > $REPO/files.sha256
echo "[$(date -Is)] done"
