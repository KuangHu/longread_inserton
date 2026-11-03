#!/bin/bash
# Submit one read-site search batch as a SLURM dependency chain (no downloads; reads from the ready manifest only):
#   [readmap array, positive test only] -> panel -> search array -> combine -> evidence -> structure -> sites
# Usage: read_site_search_submit.sh <repo batch dir> <scratch work dir>
# Job ids are appended to <repo batch dir>/jobs.txt. Acceptance, archive and cleanup run after (scale_sites.py).
set -euo pipefail
OUT=$(readlink -f $1); WORK=$2
HERE=$(dirname $(readlink -f $0))/..
N=$(($(wc -l < $OUT/samples.tsv) - 1))
mkdir -p $WORK
E="--export=ALL,OUT=$OUT,WORK=$WORK --parsable"
S=$HERE/read_site_search.sbatch
D=""
if [ -f $WORK/edited/manifest.tsv ]; then  # positive test: genome alignments against the edited assemblies first
  m=$(sbatch $E --array=0-$((N - 1))%25 --cpus-per-task=16 --mem=48G --time=6:00:00 $S readmap)
  D="--dependency=afterok:$m"; echo "readmap $m" >> $OUT/jobs.txt
fi
p=$(sbatch $E $D --cpus-per-task=4 --mem=32G $S panel)
s=$(sbatch $E --dependency=afterok:$p --array=0-$((N - 1))%25 --cpus-per-task=8 --mem=24G --time=4:00:00 $S search)
c=$(sbatch $E --dependency=afterany:$s --cpus-per-task=2 --mem=16G --time=1:00:00 $S combine)
e=$(sbatch $E --dependency=afterok:$c --cpus-per-task=8 --mem=32G $S evidence)
t=$(sbatch $E --dependency=afterok:$e --cpus-per-task=8 --mem=32G $S structure)
x=$(sbatch $E --dependency=afterok:$t --cpus-per-task=2 --mem=8G --time=1:00:00 $S sites)
echo "$(date -Is) commit $(git -C $HERE rev-parse --short HEAD) panel $p search $s combine $c evidence $e structure $t sites $x" | tee -a $OUT/jobs.txt
