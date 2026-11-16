#!/bin/bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /global/scratch/users/kh36969/longread_callout/env/msi
cd ~/middle_insertion_callout/longread_inserton/experiments/core_context
/usr/bin/time -v python ucd_combine.py combine --out scale_comamonadaceae_v0.1/ucd_batch_v0.2/final --work /global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.2/final --ids /global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.2/remaining/catalogue_ids.txt.gz --reuse /global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.1 --reuse-version ucd_batch_v0.1@faa63c3 --replace /global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.2/rerun --replace-ids scale_comamonadaceae_v0.1/ucd_batch_v0.2/rerun/affected_ids.txt --runs /global/scratch/users/kh36969/longread_callout/core_context/ucd_batch_v0.2/remaining
