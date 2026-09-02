#!/usr/bin/env bash
# The per-process progress table that NXF_ANSI_LOG=false suppresses. Nextflow's
# ANSI renderer has no throttle -- it is either continuous redraw (which is what
# made the slurm .out files 137 KB) or off -- so this reconstructs the same
# counts from the driver's log plus the work dirs. Rerun it whenever you want a
# fresh reading.
#
#   ./progress.sh                    # the current/last run
#   NXF_LOG=.nextflow.log.1 ./progress.sh    # an earlier run
#
# Reads .nextflow.log and work/*/*/.exitcode only: no nextflow, no JVM, no
# container, no ps. Safe on a login node and while the pipeline is running.
set -euo pipefail

log=${NXF_LOG:-.nextflow.log}
work=${NXF_WORK:-work}

[ -f "$log" ] || { echo "no $log here" >&2; exit 1; }

tmp=$(mktemp "${TMPDIR:-/tmp}/bea2progress.XXXXXX")
trap 'rm -f "$tmp"' EXIT

# hash -> exit status for every task dir that has finished
find "$work" -mindepth 3 -maxdepth 3 -name .exitcode \
     -exec grep -H '' {} + 2>/dev/null > "$tmp.exit" || true
trap 'rm -f "$tmp" "$tmp.exit"' EXIT

{ cat "$tmp.exit"; echo "---"; cat "$log"; } | awk '
/^---$/ { part = 2; next }
part != 2 {
    # work/5e/aca29b<...>/.exitcode:0
    if (match($0, /[0-9a-f]{2}\/[0-9a-f]{6}/)) {
        code = $0; sub(/.*:/, "", code)
        done_[substr($0, RSTART, RLENGTH)] = code
    }
    next
}
{
    if (!match($0, /\[[0-9a-f]{2}\/[0-9a-f]{6}\] (Submitted|Cached) process > [A-Za-z0-9_]+/)) next
    seg  = substr($0, RSTART, RLENGTH)
    hash = substr(seg, 2, 9)
    proc = seg; sub(/.*process > /, "", proc)
    procs[proc] = 1
    if (seg ~ /Cached/) { cached[proc]++; next }
    n[proc]++
    if (hash in done_) {
        if (done_[hash] + 0 == 0) ok[proc]++; else bad[proc]++
    } else run[proc]++
}
END {
    printf "%-26s %7s %7s %7s %7s %7s\n", "process", "done", "run", "fail", "cached", "total"
    for (p in procs) {
        t = n[p] + cached[p]
        printf "%-26s %7d %7d %7d %7d %7d\n", p, ok[p], run[p], bad[p], cached[p], t
        T_ok += ok[p]; T_run += run[p]; T_bad += bad[p]; T_c += cached[p]; T_t += t
    }
    printf "%-26s %7d %7d %7d %7d %7d\n", "TOTAL", T_ok, T_run, T_bad, T_c, T_t
}' > "$tmp"

head -1 "$tmp"
tail -n +2 "$tmp" | grep -v '^TOTAL' | sort
grep '^TOTAL' "$tmp"

grep -qE 'Goodbye|bea2 finished:' "$log" && echo "(run has finished)" || true
