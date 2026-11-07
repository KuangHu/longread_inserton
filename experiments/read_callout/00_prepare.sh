#!/bin/bash
# Transfer node only (lrc-xfer, outside SLURM): reads (ENA) and GCA/GCF assemblies (NCBI) of every sample in the
# table, so that the SLURM steps only check inputs and compute.
#   reads       scripts/reads_ready.py prepare (resume from .part, size / MD5 / FASTQ content checks, final file name
#               and verify.json status verified only after every check, ready manifest rewritten)
#   assemblies  NCBI datasets package of each accession, kept when already unpacked
# Usage (in tmux on lrc-xfer; rerun to resume): bash 00_prepare.sh   [SAMPLES=<tsv> for another panel]
set -euo pipefail
[ -z "${SLURM_JOB_ID:-}" ] || { echo "00_prepare.sh refused inside SLURM job $SLURM_JOB_ID: run it on lrc-xfer" >&2; exit 2; }
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/read_callout
DATASETS=/global/home/users/kh36969/.conda/envs/claude-env/bin/datasets
SAMPLES=${SAMPLES:-$HERE/samples.tsv}
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
fail=0
python $CODE/scripts/reads_ready.py prepare --runs $(awk -F'\t' 'NR>1{print $2}' $SAMPLES) || fail=1
while IFS=$'\t' read -r SAMPLE RUN GCA GCF PLATFORM PRESET SPLIT; do
    ASM_DIR=$PROJ/data/assemblies/$SAMPLE
    mkdir -p $ASM_DIR
    $DATASETS --version > $ASM_DIR/datasets_version.txt
    for acc in $GCA $GCF; do
        [ -d $ASM_DIR/$acc ] && continue
        # unpack into a temporary directory; the final name appears only after a complete unzip
        if $DATASETS download genome accession $acc --include genome,gbff --filename $ASM_DIR/$acc.zip --no-progressbar \
            && unzip -q -o $ASM_DIR/$acc.zip -d $ASM_DIR/$acc.tmp; then
            mv $ASM_DIR/$acc.tmp $ASM_DIR/$acc
        else
            echo "$SAMPLE $acc download failed" >&2; fail=1
        fi
    done
done < <(tail -n +2 $SAMPLES)
exit $fail
