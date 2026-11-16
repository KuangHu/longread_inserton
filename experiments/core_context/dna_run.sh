#!/bin/bash
# Driver of one DNA-entry batch (dna_discovery.py v0.2.x + ucd_batch two checks), resumable: every stage is skipped
# when its units are done; failed / missing units are resubmitted (up to RETRIES times), then the run stops with
# the missing units named. Runs on the login node (nohup); compute only through sbatch; reads verified FNA only.
#
# usage: dna_run.sh <work> <manifest> <tag> <batches> [max_pairs_per_bucket] [max_gens]
#   work        one work dir for all batches (library shards, parts and alignment cache accumulate there, so a new
#               batch is searched against the earlier ones: new library x old parts, old library x new parts)
#   batches     manifest batch numbers of this tag, e.g. 4
#   max_pairs   enumerate budget per bucket and generation (0 = all); chunks beyond it stay deferred (resumable)
#   max_gens    enumerate/align generations of this batch (deferred chunks are resumed in each next generation)
set -uo pipefail
if [ -z "${DNA_RUN_FROZEN:-}" ]; then  # run from a frozen copy: editing this file never alters a running driver
  mkdir -p ${1:?work}/runlog
  cp "$0" $1/runlog/dna_run.${3:?tag}.$(date +%s).sh
  DNA_RUN_FROZEN=1 exec bash $(ls -t $1/runlog/dna_run.$3.*.sh | head -1) "$@"
fi
WORK=${1:?work}; MAN=${2:?manifest}; TAG=${3:?tag}; BATCHES=${4:?batches}; MAXP=${5:-4000000}; MAXG=${6:-3}
CODE=/global/home/users/kh36969/middle_insertion_callout/longread_inserton/experiments/core_context
PY=/global/scratch/users/kh36969/longread_callout/env/msi/bin/python
EXCL=n0004.lr6,n0065.lr6,n0072.lr6,n0018.lr6,n0054.lr6
RETRIES=2
BIG=--partition=lr_bigmem   # single-task, memory-heavy steps (n0272, 1.5 TB)
SEED_CHUNK=${SEED_CHUNK:-75}; PART=${PART:-25}; UCD_CHUNK=${UCD_CHUNK:-5000}; EPACK=${EPACK:-8}
# DESIGN_ONLY=1 (default): a generation enumerates only the priority chunks (context sample per library sequence);
# every other pair stays deferred and is reported as pending. DESIGN_ONLY=0: all chunks within max_pairs.
DESIGN=""; [ "${DESIGN_ONLY:-1}" = 1 ] && { DESIGN="--design-only"; MAXG=1; }  # one generation holds every priority chunk
cd $CODE
mkdir -p $WORK/runlog
exec 9>$WORK/runlog/driver.lock
flock -n 9 || { echo "another dna_run.sh holds $WORK/runlog/driver.lock"; exit 1; }
LOG=$WORK/runlog/${TAG}.log
JOBS=$WORK/runlog/${TAG}_jobs.tsv
say() { echo "[$(date -Is)] $*" | tee -a $LOG; }
die() { say "STOP: $*"; exit 1; }

wait_jobs() {  # wait until none of the given job ids is pending/running (tolerates controller outages)
  local ids=$(echo "$@" | tr ' ' ',')
  [ -z "$ids" ] && return
  while true; do
    out=$(squeue -h -j $ids -o %i 2>/dev/null); rc=$?
    if [ $rc -eq 0 ] && [ -z "$out" ]; then break; fi
    sleep 60
  done
}

submit() {  # submit <stage-label> <index list "0 1 2 ..."> <sbatch args...> -- <script and args...>
  local label=$1 idx=$2; shift 2
  local sargs=() ; while [ "$1" != "--" ]; do sargs+=("$1"); shift; done; shift
  local ids=() base slice call=call$(date +%s%N)
  for base in $(for i in $idx; do echo $(( i / 1000 * 1000 )); done | sort -n | uniq); do
    slice=$(for i in $idx; do [ $(( i / 1000 * 1000 )) -eq $base ] && echo $(( i - base )); done | paste -sd,)
    local j
    for try in 1 2 3 4 5; do
      local exp="--export=ALL,OFFSET=$base" a2=()
      for x in "${sargs[@]}"; do  # a stage's own --export (ucd_batch OUT / WORK) is extended with OFFSET
        case $x in --export=*) exp="$x,OFFSET=$base";; *) a2+=("$x");; esac; done
      j=$(sbatch --parsable $exp --account=pc_rubinlab --partition=lr6 --qos=lr_normal --exclude=$EXCL \
          --output=/global/home/users/kh36969/logs/core_context/%x_%A_%a.out --job-name=d_${TAG}_${label} \
          --array=$slice "${a2[@]}" "$@" 2>>$LOG) && break
      sleep 60
    done
    [ -z "$j" ] && die "sbatch failed for $label"
    ids+=($j); echo -e "$TAG\t$label\t$j\toffset=$base\t$(date -Is)\t$call" >> $JOBS
  done
  say "$label submitted: ${ids[*]} ($(echo $idx | wc -w) tasks)"
  wait_jobs "${ids[@]}"
  # nodes where a task failed are excluded from later submissions (scratch I/O errors, node failures)
  # (only when most tasks succeeded: if most fail, the cause is the job, not the nodes)
  local bad=$(sacct -j $(echo ${ids[*]} | tr ' ' ',') -X -n -o State,NodeList%30 2>/dev/null |
              awk '{n++} $1!="COMPLETED" && $1!="RUNNING" && $1!="PENDING" && $2!="None" {f++; b[$2]=1}
                   END{if (f*4 < n) for (x in b) print x}' | sort -u | paste -sd,)
  if [ -n "$bad" ]; then EXCL="$EXCL,$bad"; say "excluding nodes after failures: $bad"; echo "$bad" >> $WORK/runlog/bad_nodes.txt; fi
}

run_until_done() {  # run_until_done <label> <missing-fn> <sbatch args...> -- <script and args...>
  local label=$1 fn=$2; shift 2
  local miss
  for r in $(seq 0 $RETRIES); do
    miss=$($fn) || die "$label: listing the missing units failed (not treated as complete)"
    [ -z "$miss" ] && { say "$label complete"; return 0; }
    [ $r -gt 0 ] && say "$label retry $r: $(echo $miss | wc -w) units"
    if [ $r -gt 0 ] && [ -n "${RETRY_ARGS:-}" ]; then   # e.g. more memory: placed after the stage's own sbatch
      local pre=() post=() seen=0 x                      # args (the later option wins), before the "--"
      for x in "$@"; do if [ $seen = 0 ] && [ "$x" = "--" ]; then seen=1; fi
        if [ $seen = 0 ]; then pre+=("$x"); else post+=("$x"); fi; done
      submit $label "$miss" "${pre[@]}" $RETRY_ARGS "${post[@]}"
    else submit $label "$miss" "$@"; fi
  done
  miss=$($fn) || die "$label: listing the missing units failed"
  [ -z "$miss" ] || die "$label units still missing after $RETRIES retries: $(echo $miss | head -c 300)"
}

# ---------------------------------------------------------------- discovery
say "batch $TAG batches=$BATCHES manifest=$MAN max_pairs=$MAXP max_gens=$MAXG commit=$(git -C $CODE rev-parse --short HEAD)"
# (awk decides alone: `awk | grep -q` fails under pipefail when grep exits before awk has written everything)
if ! awk -F'\t' -v t=$TAG '$2==t{f=1} END{exit !f}' $WORK/chunks.tsv 2>/dev/null; then
  $PY dna_discovery.py plan --work $WORK --manifest $MAN --tag $TAG --batches $BATCHES --chunk-size $SEED_CHUNK \
    --part-size $PART >> $LOG 2>&1 || die "plan failed"
fi
NCH=$(awk -F'\t' -v t=$TAG '$2==t' $WORK/chunks.tsv | wc -l)
say "seed chunks $NCH, parts $(awk -F'\t' -v t=$TAG '$2==t' $WORK/parts.tsv | wc -l)"
seeds_missing() { awk -F'\t' -v t=$TAG '$2==t{print n++" "$1}' $WORK/chunks.tsv | while read i c; do
  [ -f $WORK/seeds/$c/done.json ] || echo $i; done; }
# a chunk whose all-vs-all search needs more than a standard node (repeat-rich large genomes) is retried on lr_bigmem
RETRY_ARGS="$BIG --mem=700G"
run_until_done seeds seeds_missing --cpus-per-task=32 --mem=90G --time=24:00:00 -- \
  dna_batch.sbatch $WORK $MAN seeds --tag $TAG
RETRY_ARGS=""

if ! grep -q "\"tag\": \"$TAG\"" $WORK/library/lib_*.json 2>/dev/null; then
  submit library "0" --cpus-per-task=16 --mem=90G --time=24:00:00 -- dna_batch.sbatch $WORK $MAN library --tag $TAG
  grep -q "\"tag\": \"$TAG\"" $WORK/library/lib_*.json || die "library failed"
fi
say "library: $(grep -l "\"tag\": \"$TAG\"" $WORK/library/lib_*.json)"

occ_missing() { local o; o=$($PY dna_discovery.py occurrences --work $WORK --manifest $MAN --list) || return 1
  echo "$o" | tr ',' ' '; }
run_until_done occurrences occ_missing --cpus-per-task=32 --mem=90G --time=24:00:00 -- \
  dna_batch.sbatch $WORK $MAN occurrences

# ---------------------------------------------------------------- enumerate / align generations
touch $WORK/gens.tsv
C=$WORK/collect_$TAG
for k in $(seq 1 $MAXG); do
  [ -f $C/collect.json ] && break   # a resumed batch: its generations are closed once collected
  G=$(wc -l < $WORK/gens.tsv)
  NG=$(awk -F'\t' -v t=$TAG '$2==t && $4!="aligned_not_checked"' $WORK/gens.tsv | wc -l)
  [ $NG -ge $MAXG ] && break   # this batch already ran MAXG generations
  if [ $NG -gt 0 ]; then  # resume only if deferred chunks are left
    PREV=$((G - 1))
    pend=$(zcat $WORK/enum/g$PREV/b*/progress.tsv.gz | awk -F'\t' '$6=="deferred"' | wc -l)
    [ $pend -eq 0 ] && break
    say "generation $G: resuming $pend deferred chunks"
  fi
  # lr6 jobs are charged a whole node: the 64 buckets run $EPACK to a node job (array task t = buckets t*EPACK ..)
  enum_missing() { for t in $(seq 0 $((64 / EPACK - 1))); do for b in $(seq $((t * EPACK)) $((t * EPACK + EPACK - 1))); do
    [ -f $WORK/enum/g$G/b$(printf %02d $b)/done.json ] || { echo $t; break; }; done; done; }
  run_until_done enumerate_g$G enum_missing --cpus-per-task=32 --mem=90G --time=24:00:00 --export=ALL,PACK=$EPACK -- \
    dna_batch.sbatch $WORK $MAN enumerate --gen $G --max-pairs $MAXP $DESIGN
  align_missing() { for s in $(seq 0 127); do [ -f $WORK/align/g$G/k$(printf %03d $s)/done.json ] || echo $s; done; }
  run_until_done align_g$G align_missing --cpus-per-task=32 --mem=90G --time=24:00:00 -- \
    dna_batch.sbatch $WORK $MAN align --gen $G
  echo -e "$G\t$TAG\t$(date -Is)" >> $WORK/gens.tsv
done
# generations marked aligned_not_checked in gens.tsv (4th column) stay on disk but are not collected / checked
GENS=$(awk -F'\t' '$4!="aligned_not_checked"{print $1}' $WORK/gens.tsv | paste -sd,)
[ -f $C/collect.json ] && GENS=$($PY -c "import json;print(json.load(open('$C/collect.json'))['gens'])")

# ---------------------------------------------------------------- collect and the two checks
if [ ! -f $C/collect.json ]; then
  submit collect "0" $BIG --cpus-per-task=4 --mem=400G --time=24:00:00 -- \
    --wrap="cd $CODE; $PY dna_discovery.py collect --work $WORK --manifest $MAN --out $C --gens $GENS"
  [ -f $C/collect.json ] || die "collect failed"
fi
U=$WORK/ucd_$TAG; UW=$WORK/ucd_work_$TAG; F=$WORK/final_$TAG
PREVF=$(tail -1 $WORK/finals.txt 2>/dev/null)
if [ ! -f $U/freeze.json ]; then
  PI=""; [ -n "$PREVF" ] && PI="--previous-index $PREVF/index.tsv.gz"
  submit freeze "0" $BIG --cpus-per-task=4 --mem=400G --time=24:00:00 -- \
    --wrap="cd $CODE; $PY ucd_batch.py freeze --out $U --work $UW --manifest $MAN --input-table $C/ucd_input.tsv.gz --kinds candidate_insertion --chunk-size $UCD_CHUNK --sort-mem 100G $PI"
  [ -f $U/freeze.json ] || die "freeze failed"
fi
grep -q '"verified"' $U/freeze.json || [ "$($PY -c "import json;print(json.load(open('$U/freeze.json'))['chunks'])")" = 0 ] ||
  die "freeze verified no source FASTA ($U/sources.tsv.gz)"
NU=$($PY -c "import json;print(json.load(open('$U/freeze.json'))['chunks'])")
ucd_missing() { for c in $(seq 0 $((NU - 1))); do [ -f $UW/chunks/c$(printf %04d $c)/done.json ] || echo $c; done; }
if [ $NU -gt 0 ]; then
  run_until_done checks ucd_missing --cpus-per-task=32 --mem=64G --time=04:00:00 \
    --export=ALL,OUT=$U,WORK=$UW -- ucd_batch.sbatch chunk --manifest $MAN
  if [ ! -f $U/merge.json ]; then
    submit merge "0" $BIG --cpus-per-task=16 --mem=400G --time=24:00:00 --export=ALL,OUT=$U,WORK=$UW -- ucd_batch.sbatch merge --manifest $MAN
  fi
  grep -q '"complete": true' $U/completeness.json || die "checks incomplete: $U/completeness.json"
fi
if [ ! -f $F/combined_summary.json ]; then   # written last by combine (a table alone may be a partial file)
  RU=""; [ -n "$PREVF" ] && RU="--reuse $PREVF --reuse-pattern all_rows_%s.tsv.gz --reuse-ids $UW/reuse_ids.txt.gz"
  RN=""; [ $NU -gt 0 ] && RN="--runs $UW"
  submit combine "0" $BIG --cpus-per-task=4 --mem=300G --time=24:00:00 -- \
    --wrap="cd $CODE; $PY ucd_combine.py combine --out $F --work $F --ids $UW/catalogue_ids.txt.gz --kinds candidate_insertion $RU $RN"
  [ -f $F/combined_summary.json ] || die "combine failed"
fi
# fingerprint index of this final (sources: every code_version label with the input tables that hold its rows)
if [ ! -f $F/index.tsv.gz ]; then
  if [ $NU -gt 0 ]; then
    CV=$(zcat $UW/rows_candidate_insertion.tsv.gz | head -2 | $PY -c "import sys;h=sys.stdin.readline().rstrip('\n').split('\t');v=sys.stdin.readline().rstrip('\n').split('\t');print(dict(zip(h,v))['code_version'])")
    echo -e "$CV\t$UW/input_pairs.tsv.gz\t$U/freeze.json" >> $WORK/index_sources.tsv
  fi
  touch $WORK/index_sources.tsv
  SRC=$($PY - <<PYX
from collections import OrderedDict
d = OrderedDict()
for l in open('$WORK/index_sources.tsv'):
    cv, inp, frz = l.rstrip('\n').split('\t')
    d.setdefault(cv, [[], frz])[0].append(inp)
print(' '.join("--source '%s|%s|%s'" % (cv, ','.join(v[0]), v[1]) for cv, v in d.items()))
PYX
)
  if [ -z "$SRC" ]; then  # nothing checked in any batch so far: an empty index
    echo -e "pair_id\tinterval_kind\tfingerprint\tcode_version" | gzip > $F/index.tsv.gz
  else
  submit index "0" $BIG --cpus-per-task=4 --mem=300G --time=24:00:00 -- \
    --wrap="cd $CODE; $PY ucd_batch.py index --out $F/index.tsv.gz.tmp --final $F --kinds candidate_insertion --manifest $MAN $SRC && mv $F/index.tsv.gz.tmp $F/index.tsv.gz"
  fi
  [ -f $F/index.tsv.gz ] || die "index failed"
fi
[ -f $F/target_700_2500_summary.json ] || $PY ucd_target.py --final $F >> $LOG 2>&1 || die "target table failed"
if [ ! -f $F/core_contexts/core_contexts_summary.json ]; then
  submit contexts "0" $BIG --cpus-per-task=2 --mem=200G --time=24:00:00 -- \
    --wrap="cd $CODE; $PY ucd_core_contexts.py --final $F --dna-work $WORK --pair-libraries $(ls $WORK/collect_*/pair_other_libraries.tsv.gz 2>/dev/null | paste -sd,) --out $F/core_contexts"
  [ -f $F/core_contexts/core_contexts_summary.json ] || die "core contexts failed"
fi
grep -qx "$F" $WORK/finals.txt 2>/dev/null || echo $F >> $WORK/finals.txt

# ---------------------------------------------------------------- batch report
$PY dna_batch_report.py --work $WORK --manifest $MAN --tag $TAG --gens $GENS --collect $C --ucd $U --final $F \
  --jobs $JOBS --out $WORK/report_$TAG.json >> $LOG 2>&1 || say "report failed (results are complete)"
say "batch $TAG done: $WORK/report_$TAG.json"
