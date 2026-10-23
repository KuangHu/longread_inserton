#!/bin/bash
#SBATCH --job-name=cs1_t0
#SBATCH --account=pc_rubinlab
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/kh36969/longread_callout/logs/%x_%j.out
# T0 on a compute node: doctor, unit tests, native demo, AR20-38 frozen-input regression.
set -euo pipefail
PROJ=/global/scratch/users/kh36969/longread_callout
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton
HERE=$CODE/experiments/callout_scale_v1
OUT=$PROJ/callout_scale_v1/t0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $PROJ/env/msi
export PATH=$PROJ/env/bin:$PATH
[ ! -e "$OUT" ] || { echo "Refusing to reuse existing $OUT" >&2; exit 2; }
mkdir -p $OUT
LOG_DIR=$OUT
source $HERE/stage.sh
git -C $CODE rev-parse HEAD > $OUT/source_commit.txt

stage doctor bash -c "multisite-insertion doctor > $OUT/doctor.json"
stage unit_tests bash -c "cd $CODE && python -m unittest discover -s tests -v > $OUT/unit_tests.log 2>&1"
stage demo_data python $CODE/examples/make_demo.py --out $OUT/demo/data
stage demo_run multisite-insertion run --manifest $OUT/demo/data/manifest.json --out $OUT/demo/run \
    --methods PILER RS --subset all --threads 2
stage demo_validate bash -c "multisite-insertion validate --run $OUT/demo/run > $OUT/demo/validate.json"
stage ar2038_run multisite-insertion run --manifest $PROJ/data/assemblies/PaAR2038/manifest.json \
    --out $OUT/ar2038/callout --methods PILER RS --subset all --threads 2
stage ar2038_validate bash -c "multisite-insertion validate --run $OUT/ar2038/callout > $OUT/ar2038/validate.json"
stage ar2038_compare python $HERE/compare_callout_runs.py --old $PROJ/runs/PaAR2038/callout \
    --new $OUT/ar2038/callout --out $OUT/ar2038/regression.json
echo "[$(date -Is)] done"
