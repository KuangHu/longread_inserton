#!/bin/bash
# Transfer node only (lrc-xfer, outside SLURM): download and verify the inputs of batch10 / batch50 sample tasks,
# so that batch_sample.sbatch only checks inputs and computes.
#   reads       scripts/reads_ready.py prepare (resume from .part, size / MD5 / FASTQ content checks, final file name
#               and verify.json status verified only after every check, ready manifest rewritten)
#   assemblies  batch10/prepare_assembly.py (NCBI datasets); kept when manifest.json already exists
# Usage (in tmux on lrc-xfer; rerun the same command to resume):
#   bash prepare_inputs.sh <batch10|batch50> [order ...]     default: every row of <batch>/samples.tsv
set -euo pipefail
[ -z "${SLURM_JOB_ID:-}" ] || { echo "prepare_inputs refused inside SLURM job $SLURM_JOB_ID: run it on lrc-xfer" >&2; exit 2; }
BATCH=${1:?batch10|batch50}; shift
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/callout_scale_v1
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
# order, sample, label (manifest-safe id), run, gca, gcf of the requested rows
rows=$(python - "$HERE/$BATCH/samples.tsv" "$@" <<'PY'
import csv, re, sys
want = set(sys.argv[2:])
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    if not want or r["order"] in want:
        print("\t".join([r["order"], r["sample_id"], re.sub(r"[^A-Za-z0-9_-]", "_", r["sample_id"]), r["run"], r["gca"], r["gcf"]]))
PY
)
fail=0
python $CODE/scripts/reads_ready.py prepare --runs $(cut -f4 <<< "$rows") || fail=1
while IFS=$'\t' read -r ORDER SAMPLE LABEL RUN GCA GCF; do
    ASM=$PROJ/data/assemblies_$BATCH/$SAMPLE
    if [ -f $ASM/manifest.json ]; then echo "$SAMPLE assembly kept: $ASM"; continue; fi
    # batch10 used the sample id as the manifest id; batch50 uses the manifest-safe label
    ID=$SAMPLE; [ $BATCH = batch50 ] && ID=$LABEL
    python $HERE/batch10/prepare_assembly.py --sample $ID --gca $GCA --gcf $GCF --out $ASM || { echo "$SAMPLE assembly failed" >&2; fail=1; }
done <<< "$rows"
exit $fail
