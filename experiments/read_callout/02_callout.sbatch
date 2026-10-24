#!/bin/bash
#SBATCH --job-name=rc02_callout
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/longread_callout/logs/%x_%A_%a.out
# Assembly-only callout: native PILER/RepeatScout discovery, all families, boundaries and flanks.
source /global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments/read_callout/common.sh
mkdir -p $OUT
/usr/bin/time -v multisite-insertion run --manifest $ASM_DIR/manifest.json --out $OUT/callout \
    --subset all --threads $SLURM_CPUS_PER_TASK 2> $OUT/callout_time.txt
multisite-insertion validate --run $OUT/callout > $OUT/callout_validate.json
echo "[$(date -Is)] done"
