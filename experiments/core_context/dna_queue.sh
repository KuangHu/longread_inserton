#!/bin/bash
# Queue of DNA-entry runs, each in its own work dir: inputs -> run -> verify -> archive -> cleanup, unattended.
# Runs on the login node (nohup); compute only through sbatch (dna_run.sh). Several entries run side by side.
#
# usage: dna_queue.sh run    <root> <queue.tsv> [parallel=2]
#        dna_queue.sh status <root> <queue.tsv>
#   root       <root>/<tag>/work is the work dir of an entry, <root>/archive/<tag> its archived results,
#              <root>/queue.log the queue log, <root>/<tag>/state the last finished stage of the entry
#   queue.tsv  tab-separated, no header, '#' lines skipped:  tag  batches  manifest
#              (manifest batch numbers of the entry, e.g. 1,2,3; the manifest must list verified local FASTA files)
#
# Stages of an entry (each skipped when its state is already reached; rerunning the queue resumes):
#   inputs   every FASTA of the entry's manifest rows exists and is non-empty (content hashes are checked by the run
#            itself: seeds / occurrences and freeze verify fna_sha256). Nothing is downloaded here.
#   run      dna_run.sh (design-only generation, the two checks, target table, context tiers, report)
#   verify   report, completeness "complete": true, combined summary and context summary are all present
#   archive  small result files are copied to <root>/archive/<tag> with sha256 sums, and the sums are re-checked
#   cleanup  only after archive: re-derivable intermediates of THIS entry's work dir are listed, then deleted:
#              occ/  align/  ucd_work_<tag>/chunks/  and the pair lists under enum/ (progress.tsv.gz and done.json
#              stay, they are the record of the deferred chunks). Bytes freed go to archive/<tag>/cleanup.json.
#            Kept in place: seeds, library, collect_<tag>, ucd_<tag>, ucd_work_<tag> tables, final_<tag>, runlog, time.
# An entry that fails stops at its stage (state says where); the queue goes on with the next entry.
set -uo pipefail
MODE=${1:?run|status}; ROOT=${2:?root}; QUEUE=${3:?queue.tsv}; PAR=${4:-2}
HERE=$(cd $(dirname $0) && pwd)
PY=/global/scratch/users/kh36969/longread_callout/env/msi/bin/python
mkdir -p $ROOT/archive
QLOG=$ROOT/queue.log
say() { echo "[$(date -Is)] $*" | tee -a $QLOG; }
entries() { grep -v '^#' $QUEUE | awk -F'\t' 'NF>=3'; }

tiers() {  # one line of context tiers of a finished entry
  $PY - "$1" <<'PYX'
import json, sys
s = json.load(open(sys.argv[1]))
g = s['cores_ge']
f = lambda k: '/'.join(str(g[k].get(x, 0)) for x in ('ge4', 'ge6', 'ge8'))
print('cores>=2 %d | >=4/>=6/>=8 any %s | within one assembly %s | distinct assemblies %s | insufficient evidence %d'
      % (s['cores'], f('any'), f('single_assembly'), f('distinct_assemblies'),
         s['cores_by_evidence'].get('insufficient_evidence', 0)))
PYX
}

one() {  # one <tag> <batches> <manifest>
  local tag=$1 batches=$2 man=$3
  local D=$ROOT/$tag W=$ROOT/$tag/work A=$ROOT/archive/$tag F=$ROOT/$tag/work/final_$tag
  mkdir -p $W
  local st=$(cat $D/state 2>/dev/null)
  [ "$st" = cleaned ] && { say "$tag: already cleaned"; return 0; }
  # ---- inputs
  local miss
  miss=$($PY - "$man" "$batches" <<'PYX'
import csv, os, sys
b = set(sys.argv[2].split(','))
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t') if r['batch'] in b]
bad = [r['sample'] for r in rows if not (os.path.isfile(r['fasta']) and os.path.getsize(r['fasta']) > 0)]
print(len(rows), len(bad), ' '.join(bad[:5]))
PYX
) || { say "$tag: STOP inputs (manifest unreadable: $man)"; return 1; }
  read n nbad ex <<< "$miss"
  [ "$n" -gt 0 ] && [ "$nbad" -eq 0 ] || { say "$tag: STOP inputs ($nbad of $n FASTA files missing: $ex)"; return 1; }
  say "$tag: inputs ok ($n assemblies, batches $batches)"
  # ---- run
  if [ ! -f $W/report_$tag.json ]; then
    echo running > $D/state
    DESIGN_ONLY=1 UCD_CHUNK=20000 bash $HERE/dna_run.sh $W $man $tag $batches 0 1 > $D/driver.out 2>&1 ||
      { echo run_failed > $D/state; say "$tag: STOP run ($(tail -1 $W/runlog/$tag.log | cut -c1-200))"; return 1; }
  fi
  # ---- verify
  local U=$W/ucd_$tag
  [ -f $W/report_$tag.json ] && [ -f $F/combined_summary.json ] && [ -f $F/core_contexts/core_contexts_summary.json ] &&
    { [ ! -f $U/completeness.json ] || grep -q '"complete": true' $U/completeness.json; } ||
    { echo verify_failed > $D/state; say "$tag: STOP verify (a result file is missing or checks incomplete)"; return 1; }
  echo verified > $D/state
  say "$tag: verified | $(tiers $F/core_contexts/core_contexts_summary.json)"
  # ---- archive
  if [ ! -f $A/archive.json ]; then
    mkdir -p $A/core_contexts
    cp -p $W/report_$tag.json $W/plan_$tag.json $W/collect_$tag/collect.json $F/combined_summary.json \
       $F/target_700_2500_summary.json $F/strict_700_2500_candidate_insertion.tsv.gz $A/ &&
    cp -p $F/core_contexts/* $A/core_contexts/ &&
    { for f in freeze.json completeness.json merge.json; do [ ! -f $U/$f ] || cp -p $U/$f $A/ || false; done; } &&
    cp -p $W/runlog/$tag.log $W/runlog/${tag}_jobs.tsv $A/ ||
      { echo archive_failed > $D/state; say "$tag: STOP archive (copy failed)"; return 1; }
    (cd $A && find . -type f ! -name SHA256SUMS ! -name archive.json ! -name 'cleanup*' | sort | xargs sha256sum > SHA256SUMS &&
       sha256sum --quiet -c SHA256SUMS) || { echo archive_failed > $D/state; say "$tag: STOP archive (checksum)"; return 1; }
    # the copies equal the sources
    (cd $A && for f in report_$tag.json combined_summary.json strict_700_2500_candidate_insertion.tsv.gz; do
       s=$W/$f; [ -f $s ] || s=$F/$f; cmp -s $f $s || exit 1; done) ||
      { echo archive_failed > $D/state; say "$tag: STOP archive (copy differs from source)"; return 1; }
    echo "{\"tag\": \"$tag\", \"files\": $(wc -l < $A/SHA256SUMS), \"time\": \"$(date -Is)\", \"work\": \"$W\"}" > $A/archive.json
  fi
  echo archived > $D/state
  say "$tag: archived $A ($(wc -l < $A/SHA256SUMS) files)"
  # ---- cleanup (this entry's work dir only, fixed list, after archive)
  case $W in $ROOT/$tag/work) ;; *) say "$tag: STOP cleanup (unexpected work path $W)"; return 1;; esac
  local L=$A/cleanup_files.tsv.gz
  { for d in $W/occ $W/align $W/ucd_work_$tag/chunks; do [ -d $d ] && find $d -type f -printf '%s\t%p\n'; done
    [ -d $W/enum ] && find $W/enum -type f ! -name progress.tsv.gz ! -name done.json -printf '%s\t%p\n'; } | gzip > $L
  local nf=$(zcat $L | wc -l) bytes=$(zcat $L | awk -F'\t' '{s+=$1} END{printf "%.0f", s}')
  say "$tag: cleanup lists $nf files, $(awk -v b=$bytes 'BEGIN{printf "%.1f", b/2^30}') GiB ($L)"
  zcat $L | cut -f2- | while IFS= read -r f; do
    case $f in $W/occ/*|$W/align/*|$W/ucd_work_$tag/chunks/*|$W/enum/*) rm -f -- "$f";; *) echo "skip $f" >> $QLOG;; esac
  done
  for d in $W/occ $W/align $W/ucd_work_$tag/chunks; do [ -d $d ] && find $d -depth -type d -empty -delete; done
  local left=$(zcat $L | cut -f2- | while IFS= read -r f; do [ -e "$f" ] && echo x; done | wc -l)
  echo "{\"tag\": \"$tag\", \"files_listed\": $nf, \"bytes_freed\": $bytes, \"files_left\": $left, \"time\": \"$(date -Is)\"}" > $A/cleanup.json
  [ "$left" -eq 0 ] || { echo cleanup_incomplete > $D/state; say "$tag: STOP cleanup ($left listed files still there)"; return 1; }
  echo cleaned > $D/state
  say "$tag: cleaned, $(awk -v b=$bytes 'BEGIN{printf "%.1f", b/2^30}') GiB freed"
}

if [ "$MODE" = status ]; then
  entries | while IFS=$'\t' read -r tag batches man; do
    st=$(cat $ROOT/$tag/state 2>/dev/null || echo "-")
    S=$ROOT/$tag/work/final_$tag/core_contexts/core_contexts_summary.json; [ -f $S ] || S=$ROOT/archive/$tag/core_contexts/core_contexts_summary.json
    if [ -f $S ]; then echo "$tag [$st] $(tiers $S)"
    else echo "$tag [$st] $(tail -1 $ROOT/$tag/work/runlog/$tag.log 2>/dev/null | cut -c1-110)"; fi
  done
  exit 0
fi
[ "$MODE" = run ] || { echo "mode: run | status"; exit 1; }
exec 8>$ROOT/$(basename $QUEUE).lock   # one runner per queue file; entries of different files must differ
flock -n 8 || { echo "another dna_queue.sh run holds the lock of $QUEUE"; exit 1; }
say "queue start: $(entries | wc -l) entries, $PAR at a time, watch=${WATCH:-0}, code $(git -C $HERE rev-parse --short HEAD)"
# WATCH=1: the queue file is read again every 10 minutes, so entries appended later (e.g. after a download was
# verified) are picked up; touch <root>/queue.stop to let it end after the running entries.
declare -A started
while true; do
  while IFS=$'\t' read -r tag batches man; do
    [ -n "${started[$tag]:-}" ] && continue
    while [ $(jobs -rp | wc -l) -ge $PAR ]; do wait -n; done
    started[$tag]=1
    one $tag $batches $man 8>&- < /dev/null &
    sleep 5
  done < <(entries)
  [ "${WATCH:-0}" = 1 ] && [ ! -f $ROOT/queue.stop ] || break
  sleep 600
done
wait
say "queue end: $(entries | while IFS=$'\t' read -r tag b m; do echo -n "$tag=$(cat $ROOT/$tag/state 2>/dev/null) "; done)"
