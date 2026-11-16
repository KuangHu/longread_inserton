#!/bin/bash
#SBATCH --job-name=cs1_filled
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-5%2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/global/home/users/kh36969/logs/callout_scale_v1/%x_%A_%a.out
# Filled-site read check only (protocol 3.5 / 4): explicit inputs from execution_manifest.json,
# new output directory, audited classifier and checkers. No non-reference insertion step,
# no deletion of earlier evidence.
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/callout_scale_v1
# Some lr6 nodes had a hung /global/scratch (writes blocked, nothing visible elsewhere);
# probe it with a timeout before doing anything, and leave the log on /global/home.
echo "[$(date -Is)] host=$(hostname) task=$SLURM_ARRAY_TASK_ID"
if ! timeout 60 bash -c "echo ok > $PROJ/logs/.probe_${SLURM_JOB_ID} && cat $PROJ/logs/.probe_${SLURM_JOB_ID} > /dev/null"; then
    echo "scratch_unavailable on $(hostname); exiting 75 for resubmission elsewhere" >&2
    exit 75
fi
rm -f $PROJ/logs/.probe_${SLURM_JOB_ID}
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi

eval "$(python - "$HERE/execution_manifest.json" "$SLURM_ARRAY_TASK_ID" <<'EOF'
import json, shlex, sys
e = [s for s in json.load(open(sys.argv[1]))["samples"] if s["order"] == int(sys.argv[2])][0]
for k in ("sample_id", "run", "fasta", "fasta_sha256", "reads", "reads_sha256", "preset",
          "v1_sample_dir", "v1_callout", "v1_evidence", "independence", "assembly_consistency"):
    print("%s=%s" % (k.upper(), shlex.quote(str(e[k]))))
EOF
)"
OUT=$PROJ/callout_scale_v1/$SAMPLE_ID
[ ! -e "$OUT" ] || { echo "Refusing to reuse existing $OUT" >&2; exit 2; }
mkdir -p $OUT/mapping $OUT/qc
LOG_DIR=$OUT
source $HERE/stage.sh
echo "[$(date -Is)] sample=$SAMPLE_ID run=$RUN host=$(hostname)"
git -C $CODE rev-parse HEAD > $OUT/source_commit.txt
minimap2 --version > $OUT/mapping/minimap2_version.txt
MEMBERS=$V1_CALLOUT/results/member_boundaries.json

check_hash() { [ "$(sha256sum "$1" | cut -d' ' -f1)" = "$2" ] || { echo "SHA256 mismatch: $1" >&2; return 1; }; }
stage verify_inputs bash -c "$(declare -f check_hash); check_hash $READS $READS_SHA256 && check_hash $FASTA $FASTA_SHA256"

MM="-x $PRESET -c --secondary=yes -N 50 -p 0.8 -I 8G"
echo "minimap2 $MM -t $SLURM_CPUS_PER_TASK $FASTA $READS" > $OUT/mapping/command.txt
stage map bash -c "minimap2 $MM -t $SLURM_CPUS_PER_TASK $FASTA $READS 2> $OUT/mapping/minimap2.log | gzip -1 > $OUT/mapping/reads.paf.gz"

stage cap_reads python $CODE/scripts/audit_filled_site.py cap-reads --members $MEMBERS \
    --paf $OUT/mapping/reads.paf.gz --out $OUT/mapping/cap_reads.txt
RECHECK=()
if [ -s $OUT/mapping/cap_reads.txt ]; then
    stage cap_extract python $CODE/scripts/audit_filled_site.py extract --reads $OUT/mapping/cap_reads.txt \
        --fastq $READS --out $OUT/mapping/cap_reads.fastq
    MM2="-x $PRESET -c --secondary=yes -N 200 -p 0.5 -I 8G"
    echo "minimap2 $MM2 -t $SLURM_CPUS_PER_TASK $FASTA cap_reads.fastq" > $OUT/mapping/recheck_command.txt
    stage recheck_map bash -c "minimap2 $MM2 -t $SLURM_CPUS_PER_TASK $FASTA $OUT/mapping/cap_reads.fastq 2> $OUT/mapping/minimap2_recheck.log > $OUT/mapping/recheck.paf"
    RECHECK=(--recheck-paf $OUT/mapping/recheck.paf)
fi

stage classify python $CODE/scripts/audit_filled_site.py classify --members $MEMBERS --reference $FASTA \
    --paf $OUT/mapping/reads.paf.gz "${RECHECK[@]}" --v1 $V1_EVIDENCE/member_read_support.tsv \
    --sample $SAMPLE_ID --independence $INDEPENDENCE --assembly-consistency $ASSEMBLY_CONSISTENCY \
    --out $OUT/evidence_audited
stage audited_check python $CODE/scripts/check_audited_evidence.py --audit-dir $OUT/evidence_audited \
    --members $MEMBERS --reference $FASTA --out $OUT/qc/audited_check.json
stage bundle_v1_live python $CODE/scripts/check_callout_bundle.py --sample-dir $V1_SAMPLE_DIR \
    --reads-verify $PROJ/data/reads/$RUN/verify.json --out $OUT/qc/bundle_v1.json
echo "[$(date -Is)] done"
