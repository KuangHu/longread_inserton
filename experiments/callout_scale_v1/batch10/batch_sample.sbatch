#!/bin/bash
#SBATCH --job-name=cs1_batch10
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-10%2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/global/home/users/kh36969/logs/callout_scale_v1/%x_%A_%a.out
# Controlled 10-sample batch (protocol T2): one sample per task, at most 2 concurrent.
# input check (reads + assembly prepared on lrc-xfer by prepare_inputs.sh) -> GCA/GCF assembly choice -> callout (--methods PILER RS --subset all)
# -> validate -> map (fixed ONT command) -> audited filled-site check (nearby_competition_v2) -> checker.
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
export PATH=$PROJ/env/bin:$PATH

row=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" '$1==n' $HERE/batch10/samples.tsv)
IFS=$'\t' read -r ORDER SAMPLE RUN BIOSAMPLE GCA GCF PRESET <<< "$row"
# Input check before any output exists, so a task stopped with input_not_ready (exit 3) can simply be resubmitted.
python $CODE/scripts/reads_ready.py check --runs $RUN > /dev/null || exit 3
[ -f $PROJ/data/assemblies_batch10/$SAMPLE/manifest.json ] || { echo "input_not_ready assembly $SAMPLE not prepared" >&2; exit 3; }
OUT=$PROJ/batch10/$SAMPLE
[ ! -e "$OUT" ] || { echo "Refusing to reuse existing $OUT" >&2; exit 2; }
mkdir -p $OUT/mapping $OUT/qc
LOG_DIR=$OUT
source $HERE/stage.sh
git -C $CODE rev-parse HEAD > $OUT/source_commit.txt
echo "sample=$SAMPLE run=$RUN gca=$GCA gcf=$GCF"

# Inputs are prepared on the transfer node (prepare_inputs.sh); this job never downloads. A run that is not in the
# ready manifest as verified, or whose file is missing or changed, stops the job with input_not_ready (exit 3).
stage reads_ready python $CODE/scripts/reads_ready.py check --runs $RUN
READS=$(python $CODE/scripts/reads_ready.py check --runs $RUN | cut -f3 -d' ')
ASM=$PROJ/data/assemblies_batch10/$SAMPLE
# Prepared with prepare_inputs.sh; reuse only if the chosen FASTA still has the recorded hash.
[ -f $ASM/manifest.json ] || { echo "input_not_ready assembly $ASM/manifest.json missing: run prepare_inputs.sh batch10 on lrc-xfer" >&2; exit 3; }
stage assembly_ready python -c "import hashlib,json,sys; c=json.load(open('$ASM/assembly_check.json')); sys.exit(3 * (hashlib.sha256(open(c['fasta'],'rb').read()).hexdigest()!=c['fasta_sha256']))"
FASTA=$(python -c "import json;print(json.load(open('$ASM/manifest.json'))['samples'][0]['fasta'])")

# Discovery budget per protocol: BLASTX threads 2; native tools use their own defaults.
stage warm_tools bash -c "for t in blastn blastx makeblastdb; do \$t -version > /dev/null 2>&1 || true; done"
stage callout multisite-insertion run --manifest $ASM/manifest.json --out $OUT/callout --methods PILER RS --subset all --threads 2
stage validate bash -c "multisite-insertion validate --run $OUT/callout > $OUT/qc/validate.json"

MM="-x $PRESET -c --secondary=yes -N 50 -p 0.8 -I 8G"
minimap2 --version > $OUT/mapping/minimap2_version.txt
echo "minimap2 $MM -t $SLURM_CPUS_PER_TASK $FASTA $READS" > $OUT/mapping/command.txt
stage map bash -c "minimap2 $MM -t $SLURM_CPUS_PER_TASK $FASTA $READS 2> $OUT/mapping/minimap2.log | gzip -1 > $OUT/mapping/reads.paf.gz"
MEMBERS=$OUT/callout/results/member_boundaries.json
stage cap_reads python $CODE/scripts/audit_filled_site.py cap-reads --members $MEMBERS --paf $OUT/mapping/reads.paf.gz --out $OUT/mapping/cap_reads.txt
RECHECK=()
if [ -s $OUT/mapping/cap_reads.txt ]; then
    stage cap_extract python $CODE/scripts/audit_filled_site.py extract --reads $OUT/mapping/cap_reads.txt --fastq $READS --out $OUT/mapping/cap_reads.fastq
    MM2="-x $PRESET -c --secondary=yes -N 200 -p 0.5 -I 8G"
    echo "minimap2 $MM2 -t $SLURM_CPUS_PER_TASK $FASTA cap_reads.fastq" > $OUT/mapping/recheck_command.txt
    stage recheck_map bash -c "minimap2 $MM2 -t $SLURM_CPUS_PER_TASK $FASTA $OUT/mapping/cap_reads.fastq 2> $OUT/mapping/minimap2_recheck.log > $OUT/mapping/recheck.paf"
    RECHECK=(--recheck-paf $OUT/mapping/recheck.paf)
fi
stage classify python $CODE/scripts/audit_filled_site.py classify --members $MEMBERS --reference $FASTA \
    --paf $OUT/mapping/reads.paf.gz "${RECHECK[@]}" --sample $SAMPLE \
    --independence independence_unverified --assembly-consistency see_assembly_method --out $OUT/evidence_audited
stage audited_check python $CODE/scripts/check_audited_evidence.py --audit-dir $OUT/evidence_audited \
    --members $MEMBERS --reference $FASTA --out $OUT/qc/audited_check.json
du -sb $OUT $PROJ/data/reads/$RUN $ASM > $OUT/qc/disk_usage.txt
echo "[$(date -Is)] done"
