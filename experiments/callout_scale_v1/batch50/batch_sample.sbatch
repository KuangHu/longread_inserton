#!/bin/bash
#SBATCH --job-name=cs1_batch50
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/logs/callout_scale_v1/%x_%A_%a.out
# Frozen 50-sample panel, one sample per array task (task id = samples.tsv order).
# Submit in groups of ten, at most 2 concurrent, e.g.
#   sbatch --array=1-10%2 --exclude=<bad nodes> batch_sample.sbatch
# fetch/verify reads -> GCA/GCF choice + expected SHA256 -> callout (--methods PILER RS --subset all)
# -> validate -> [fixed-seed 80x downsample if read_policy=downsample] -> map (fixed ONT command)
# -> audited filled-site check (nearby_competition_v2) -> checker.
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/callout_scale_v1
echo "[$(date -Is)] host=$(hostname) task=$SLURM_ARRAY_TASK_ID"
if ! timeout 60 bash -c "echo ok > $PROJ/logs/.probe_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} && cat $PROJ/logs/.probe_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} > /dev/null"; then
    echo "scratch_unavailable on $(hostname); exiting 75 for resubmission elsewhere" >&2
    exit 75
fi
rm -f $PROJ/logs/.probe_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
export PATH=$PROJ/env/bin:$PATH

# Read the row by column name; the table has more columns than batch10.
col() { python -c "import csv,sys; r=[x for x in csv.DictReader(open('$HERE/batch50/samples.tsv'),delimiter='\t') if x['order']=='$SLURM_ARRAY_TASK_ID']; assert len(r)==1; print(r[0]['$1'])"; }
SAMPLE=$(col sample_id); RUN=$(col run); GCA=$(col gca); GCF=$(col gcf); PRESET=$(col preset)
EXPECT_SHA=$(col expected_fasta_sha256); EXPECT_SRC=$(col selected_fasta_source); POLICY=$(col read_policy)
ASM_BP=$(col assembly_bp); DS_DEPTH=$(col downsample_depth); DS_SEED=$(col downsample_seed)
OUT=$PROJ/batch50/$SAMPLE
# Package manifest IDs allow only [A-Za-z0-9_-]; the panel sample_id is kept for directories and tables.
LABEL=$(echo "$SAMPLE" | sed 's/[^A-Za-z0-9_-]/_/g')
[ ! -e "$OUT" ] || { echo "Refusing to reuse existing $OUT" >&2; exit 2; }
mkdir -p $OUT/mapping $OUT/qc
LOG_DIR=$OUT
source $HERE/stage.sh
git -C $CODE rev-parse HEAD > $OUT/source_commit.txt
echo "sample=$SAMPLE label=$LABEL run=$RUN gca=$GCA gcf=$GCF policy=$POLICY"

if [ -f $PROJ/data/reads/$RUN/verify.json ]; then
    stage reads_verified python -c "import json,sys; d=json.load(open('$PROJ/data/reads/$RUN/verify.json')); sys.exit(d['status']!='verified')"
else
    stage fetch_reads python $CODE/scripts/fetch_long_reads.py fetch --out $PROJ/data/reads --runs $RUN
fi
READS=$(python -c "import json;print(json.load(open('$PROJ/data/reads/$RUN/verify.json'))['local_file'])")
ASM=$PROJ/data/assemblies_batch50/$SAMPLE
if [ -f $ASM/manifest.json ]; then
    stage assembly_reused python -c "import hashlib,json,sys; c=json.load(open('$ASM/assembly_check.json')); sys.exit(hashlib.sha256(open(c['fasta'],'rb').read()).hexdigest()!=c['fasta_sha256'])"
else
    stage prepare_assembly python $HERE/batch10/prepare_assembly.py --sample $LABEL --gca $GCA --gcf $GCF --out $ASM
fi
# The FASTA chosen on the cluster must be the one frozen in samples.tsv.
stage assembly_expected python -c "import json,sys; c=json.load(open('$ASM/assembly_check.json')); src='gcf_sequence_and_annotation' if '$EXPECT_SRC'.startswith('gcf') else 'gca_sequence'; ok=c['fasta_sha256']=='$EXPECT_SHA' and c['choice'].startswith(src[:3]); print(c['choice'], c['fasta_sha256']); sys.exit(not ok)"
FASTA=$(python -c "import json;print(json.load(open('$ASM/manifest.json'))['samples'][0]['fasta'])")

# Discovery budget per protocol: BLASTX threads 2; native tools use their own defaults.
stage warm_tools bash -c "for t in blastn blastx makeblastdb; do \$t -version > /dev/null 2>&1 || true; done"
stage callout multisite-insertion run --manifest $ASM/manifest.json --out $OUT/callout --methods PILER RS --subset all --threads 2
stage validate bash -c "multisite-insertion validate --run $OUT/callout > $OUT/qc/validate.json"

# Reads for mapping: the full verified FASTQ, or a fixed-seed per-molecule downsample of it.
MAP_READS=$READS
if [ "$POLICY" = downsample ]; then
    DS=$PROJ/data/reads_downsampled/$RUN/${RUN}_${DS_DEPTH}x_seed${DS_SEED}
    mkdir -p $(dirname $DS)
    if [ -f $DS.json ]; then
        stage downsample_reused python -c "import hashlib,json,sys; d=json.load(open('$DS.json')); h=hashlib.sha256(open(d['output_file'],'rb').read()).hexdigest(); sys.exit(h!=d['output_sha256'])"
    else
        stage downsample python $HERE/batch50/downsample_reads.py --verify $PROJ/data/reads/$RUN/verify.json \
            --assembly-bp $ASM_BP --depth $DS_DEPTH --seed $DS_SEED --out $DS
    fi
    cp $DS.json $OUT/mapping/downsample.json
    MAP_READS=$DS.fastq.gz
fi

MM="-x $PRESET -c --secondary=yes -N 50 -p 0.8 -I 8G"
minimap2 --version > $OUT/mapping/minimap2_version.txt
echo "minimap2 $MM -t $SLURM_CPUS_PER_TASK $FASTA $MAP_READS" > $OUT/mapping/command.txt
stage map bash -c "minimap2 $MM -t $SLURM_CPUS_PER_TASK $FASTA $MAP_READS 2> $OUT/mapping/minimap2.log | gzip -1 > $OUT/mapping/reads.paf.gz"
MEMBERS=$OUT/callout/results/member_boundaries.json
stage cap_reads python $CODE/scripts/audit_filled_site.py cap-reads --members $MEMBERS --paf $OUT/mapping/reads.paf.gz --out $OUT/mapping/cap_reads.txt
RECHECK=()
if [ -s $OUT/mapping/cap_reads.txt ]; then
    stage cap_extract python $CODE/scripts/audit_filled_site.py extract --reads $OUT/mapping/cap_reads.txt --fastq $MAP_READS --out $OUT/mapping/cap_reads.fastq
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
