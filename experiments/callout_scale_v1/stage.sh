# Stage runner for callout_scale_v1 jobs; sourced after `set -euo pipefail`.
# stage NAME CMD...  runs CMD under /usr/bin/time -v, appends one row to $LOG_DIR/stages.tsv
# (name, start, end, exit code, wall, peak RSS kB, command) and stops the job on failure.
stage() {
    local name=$1; shift
    local start rc end
    start=$(date -Is)
    if /usr/bin/time -v -o "$LOG_DIR/time_$name.txt" "$@"; then rc=0; else rc=$?; fi
    end=$(date -Is)
    local wall rss
    wall=$(grep -m1 "Elapsed (wall clock)" "$LOG_DIR/time_$name.txt" | awk '{print $NF}' || true)
    rss=$(grep -m1 "Maximum resident" "$LOG_DIR/time_$name.txt" | awk '{print $NF}' || true)
    [ -s "$LOG_DIR/stages.tsv" ] || printf 'stage\tstart\tend\texit_code\twall\tmax_rss_kb\tcommand\n' > "$LOG_DIR/stages.tsv"
    local cmd
    cmd=$(printf '%s ' "$@" | tr '\n\t' '  ' | tr -s ' ')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$start" "$end" "$rc" "$wall" "$rss" "$cmd" >> "$LOG_DIR/stages.tsv"
    if [ "$rc" -ne 0 ]; then
        echo "stage $name failed with exit code $rc" >&2
        exit "$rc"
    fi
}
