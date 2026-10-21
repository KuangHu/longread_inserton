#!/usr/bin/env bash
set -euo pipefail
: "${PILOT_MINIMAP2:?Set PILOT_MINIMAP2 to the actual minimap2 binary}"
# map-hifi does not change mini_batch_size. -K affects batching, not thresholds.
exec "$PILOT_MINIMAP2" -K 50M "$@"
