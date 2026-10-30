#!/bin/bash
#SBATCH --job-name=rc03_reads
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/kh36969/longread_callout/logs/%x_%A_%a.out
# Read evidence: molecule-named reads -> minimap2 to the complete assembly -> classification.
source /global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments/read_callout/common.sh
R=$OUT/reads
mkdir -p $R
FASTQ=$(python $CODE/scripts/reads_ready.py check --runs $RUN | cut -f3 -d' ') || exit 3  # ready manifest only; never downloads
REF=$(python -c "import json;print(json.load(open('$ASM_DIR/manifest.json'))['samples'][0]['fasta'])")
minimap2 --version > $R/minimap2_version.txt

if [ ! -s $R/molecules.fastq.gz ]; then
    split=""; [ "$SPLIT" = yes ] && split="--split-adapters"
    python $HERE/molecule_reads.py rename $FASTQ $split --stats $R/rename_stats.json | gzip -1 > $R/molecules.fastq.gz.tmp
    mv $R/molecules.fastq.gz.tmp $R/molecules.fastq.gz
fi

# Secondary alignments kept (C1); -N 50 above the default of 5; single-part index for a bacterial genome.
MM="-x $PRESET -c --secondary=yes -N 50 -p 0.8 -I 8G"
echo "minimap2 $MM" > $R/minimap2_command.txt
if [ ! -s $R/reads.paf.gz ]; then
    /usr/bin/time -v minimap2 $MM -t $SLURM_CPUS_PER_TASK $REF $R/molecules.fastq.gz 2> $R/minimap2.log | gzip -1 > $R/reads.paf.gz.tmp
    mv $R/reads.paf.gz.tmp $R/reads.paf.gz
fi

python $HERE/classify_reads.py insertions --paf $R/reads.paf.gz --out $R/insertions.tsv
python $HERE/classify_reads.py extract --insertions $R/insertions.tsv --fastq $R/molecules.fastq.gz --out $R/insertions.fa
minimap2 $MM -t $SLURM_CPUS_PER_TASK $REF $R/insertions.fa 2> $R/minimap2_insertions.log > $R/insertions.paf

rm -rf $OUT/evidence
python $HERE/classify_reads.py classify --run $OUT/callout --reference $REF --paf $R/reads.paf.gz \
    --insertions $R/insertions.tsv --ins-paf $R/insertions.paf --out $OUT/evidence --sample $SAMPLE
echo "[$(date -Is)] done"
