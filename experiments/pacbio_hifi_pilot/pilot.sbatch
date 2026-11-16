#!/bin/bash
#SBATCH --job-name=hifi_pilot
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/pacbio_hifi_pilot/%x_%j.out
# PacBio HiFi single-sample compatibility pilot (Sam8H1, SRR13837247, GCA_016904925.1).
# Same route as callout_scale_v1/batch50 except: minimap2 -x map-hifi for the main mapping and the
# cap recheck; read QV statistics and an assembly/report check added; full reads, no downsampling.
# Inputs first, on lrc-xfer: bash experiments/pacbio_hifi_pilot/prepare_inputs.sh
# Submit: sbatch --exclude=<bad nodes> experiments/pacbio_hifi_pilot/pilot.sbatch
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/pacbio_hifi_pilot
SAMPLE=Sam8H1; RUN=SRR13837247; GCA=GCA_016904925.1; GCF=GCF_016904925.1
BIOSAMPLE=SAMN17825431; TOTAL_BP=4131713; PRESET=map-hifi
echo "[$(date -Is)] host=$(hostname) job=$SLURM_JOB_ID"
if ! timeout 60 bash -c "echo ok > $PROJ/logs/.probe_${SLURM_JOB_ID} && cat $PROJ/logs/.probe_${SLURM_JOB_ID} > /dev/null"; then
    echo "scratch_unavailable on $(hostname); exiting 75 for resubmission elsewhere" >&2
    exit 75
fi
rm -f $PROJ/logs/.probe_${SLURM_JOB_ID}
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
export PATH=$PROJ/env/bin:$PATH

# Input check before any output exists, so a task stopped with input_not_ready (exit 3) can simply be resubmitted.
python $CODE/scripts/reads_ready.py check --runs $RUN > /dev/null || exit 3
[ -f $PROJ/data/assemblies_pacbio_pilot/$SAMPLE/manifest.json ] || { echo "input_not_ready assembly $SAMPLE not prepared" >&2; exit 3; }
OUT=$PROJ/pacbio_hifi_pilot/$SAMPLE
[ ! -e "$OUT" ] || { echo "Refusing to reuse existing $OUT" >&2; exit 2; }
mkdir -p $OUT/mapping $OUT/qc
LOG_DIR=$OUT
source $CODE/experiments/callout_scale_v1/stage.sh
git -C $CODE rev-parse HEAD > $OUT/source_commit.txt
{ minimap2 --version | sed 's/^/minimap2 /'; multisite-insertion --version 2>&1 || true; python --version
  pals 2>&1 | head -2 || true; piler2 --version 2>&1 | head -1 || true; RepeatScout 2>&1 | head -1 || true
  blastn -version | head -1; ~/.conda/envs/claude-env/bin/datasets --version; } > $OUT/qc/tool_versions.txt 2>&1 || true

# Inputs are prepared on the transfer node (prepare_inputs.sh); this job never downloads. A run that is not in the
# ready manifest as verified, or whose file is missing or changed, stops the job with input_not_ready (exit 3).
stage reads_ready python $CODE/scripts/reads_ready.py check --runs $RUN
cp $PROJ/data/reads/$RUN/verify.json $OUT/qc/reads_verify.json
READS=$(python $CODE/scripts/reads_ready.py check --runs $RUN | cut -f3 -d' ')
stage read_stats python $HERE/read_stats.py --verify $PROJ/data/reads/$RUN/verify.json --out $OUT/qc/read_stats.json

ASM=$PROJ/data/assemblies_pacbio_pilot/$SAMPLE
# Prepared and verified on lrc-xfer (prepare_inputs.sh): only the recorded hash is checked here.
[ -f $ASM/manifest.json ] && [ -f $ASM/assembly_verify.json ] || { echo "input_not_ready assembly $ASM not prepared: run prepare_inputs.sh on lrc-xfer" >&2; exit 3; }
stage assembly_ready python -c "import hashlib,json,sys; c=json.load(open('$ASM/assembly_check.json')); sys.exit(3 * (hashlib.sha256(open(c['fasta'],'rb').read()).hexdigest()!=c['fasta_sha256']))"
cp $ASM/assembly_check.json $ASM/assembly_verify.json $OUT/qc/
FASTA=$(python -c "import json;print(json.load(open('$ASM/manifest.json'))['samples'][0]['fasta'])")

# Discovery budget per protocol: BLASTX threads 2; native tools use their own defaults. Reads are not used here.
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
# The assembly was built from these same reads: this is a compatibility/consistency check, not independent validation.
stage classify python $CODE/scripts/audit_filled_site.py classify --members $MEMBERS --reference $FASTA \
    --paf $OUT/mapping/reads.paf.gz "${RECHECK[@]}" --sample $SAMPLE \
    --independence independence_unverified --assembly-consistency assembly_built_from_same_reads --out $OUT/evidence_audited
stage audited_check python $CODE/scripts/check_audited_evidence.py --audit-dir $OUT/evidence_audited \
    --members $MEMBERS --reference $FASTA --out $OUT/qc/audited_check.json
du -sb $OUT $PROJ/data/reads/$RUN $ASM > $OUT/qc/disk_usage.txt
echo "[$(date -Is)] done"
