#!/usr/bin/env bash
# Run from a workspace containing the data directory and installed tools.
set -euo pipefail
ROOT="${1:?Usage: bash run.sh /absolute/workspace}"
export PILOT_MINIMAP2="$ROOT/minimap2_build/minimap2"
DATA="$ROOT/hifi_marine_pilot"
READS="$DATA/reads/SRR28395891_subreads.fastq.gz"
# Caller parameters are unchanged from the prototype defaults.
"$ROOT/site_recovery_venv/bin/insertion-site-recovery" run \
  --reference "$DATA/reference/marine_mag_pool.fasta" \
  --reads "$READS" --platform hifi --sample SRR28395891 \
  --threads 4 --minimap2 "$ROOT/longread_inserton/validation/real_hifi_pilot/minimap2_lowmem.sh" \
  --sniffles "$ROOT/site_recovery_venv/bin/sniffles" \
  --out "$DATA/calls"
