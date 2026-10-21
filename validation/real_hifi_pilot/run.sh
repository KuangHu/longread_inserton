#!/usr/bin/env bash
# Reproduce the final staged workflow after fetching reference and complete reads.
set -euo pipefail
ROOT="${1:?Usage: bash run.sh /absolute/workspace}"
PY="$ROOT/site_recovery_venv/bin/python"
TOOLS="$ROOT/longread_inserton/validation/real_hifi_pilot"
DATA="$ROOT/hifi_marine_pilot"
"$PY" "$TOOLS/align.py" "$DATA" "$ROOT/minimap2_build/minimap2"
"$PY" "$TOOLS/summarize.py" "$DATA" --coverage-only
"$PY" "$TOOLS/call_targets.py" "$DATA" "$ROOT/site_recovery_venv/bin/sniffles" "$ROOT/site_recovery_venv/bin/insertion-site-recovery"
"$PY" "$TOOLS/summarize.py" "$DATA"
