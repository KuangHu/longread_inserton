#!/bin/bash
#SBATCH --job-name=cs1_reclass
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-5%2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/global/home/users/kh36969/logs/callout_scale_v1/%x_%A_%a.out
# Reclassify the five samples with classifier_revision nearby_competition_v2 (NEARBY_COMPETITION_REVIEW.md):
# same PAF (hash-checked), same reference and member_boundaries; new output root; old audit kept.
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/callout_scale_v1
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
for k in ("sample_id", "run", "fasta", "fasta_sha256", "v1_sample_dir", "v1_callout", "v1_evidence",
          "independence", "assembly_consistency"):
    print("%s=%s" % (k.upper(), shlex.quote(str(e[k]))))
EOF
)"
OLD=$PROJ/callout_scale_v1/$SAMPLE_ID
OUT=$PROJ/callout_scale_v1_v2/$SAMPLE_ID
[ ! -e "$OUT" ] || { echo "Refusing to reuse existing $OUT" >&2; exit 2; }
mkdir -p $OUT/qc $OUT/mapping
LOG_DIR=$OUT
source $HERE/stage.sh
git -C $CODE rev-parse HEAD > $OUT/source_commit.txt
PAF=$OLD/mapping/reads.paf.gz
PAF_SHA=$(awk -F'\t' -v s=$SAMPLE_ID '$1==s{print $2}' $HERE/v2_input_paf_hashes.tsv)
MEMBERS=$V1_CALLOUT/results/member_boundaries.json

check_hash() { [ "$(sha256sum "$1" | cut -d' ' -f1)" = "$2" ] || { echo "SHA256 mismatch: $1" >&2; return 1; }; }
stage verify_inputs bash -c "$(declare -f check_hash); check_hash $PAF $PAF_SHA && check_hash $FASTA $FASTA_SHA256"
stage cap_reads python $CODE/scripts/audit_filled_site.py cap-reads --members $MEMBERS --paf $PAF --out $OUT/mapping/cap_reads.txt
RECHECK=()
if [ -s $OUT/mapping/cap_reads.txt ]; then
    echo "cap reads present: the recheck path from filled_site.sbatch is required; stopping" >&2
    exit 3
fi
stage classify python $CODE/scripts/audit_filled_site.py classify --members $MEMBERS --reference $FASTA \
    --paf $PAF "${RECHECK[@]}" --v1 $V1_EVIDENCE/member_read_support.tsv \
    --sample $SAMPLE_ID --independence $INDEPENDENCE --assembly-consistency $ASSEMBLY_CONSISTENCY \
    --out $OUT/evidence_audited
stage audited_check python $CODE/scripts/check_audited_evidence.py --audit-dir $OUT/evidence_audited \
    --members $MEMBERS --reference $FASTA --out $OUT/qc/audited_check.json
stage bundle_v1_live python $CODE/scripts/check_callout_bundle.py --sample-dir $V1_SAMPLE_DIR \
    --reads-verify $PROJ/data/reads/$RUN/verify.json --out $OUT/qc/bundle_v1.json
stage diff_vs_previous python $HERE/diff_audits.py --old $OLD/evidence_audited --new $OUT/evidence_audited \
    --out $OUT/qc/diff_vs_previous_audit.tsv
echo "[$(date -Is)] done"
