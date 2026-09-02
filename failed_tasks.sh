#!/usr/bin/env bash
# Per-process status tally and failure details for a run. Two things it stands
# in for:
#
#   * the "[100%] n of m" table, which NXF_ANSI_LOG=false suppresses (that flag
#     is set in run_test.sh because in a non-TTY the ANSI renderer writes every
#     redraw frame -- it is what made the old slurm .out files 137 KB). Unset it
#     to get the live table back at the cost of log size.
#   * -with-trace, which cannot be used here at all: it injects a ps-based
#     collector into each container and none of our images ship procps, so every
#     containerised task would exit 1 immediately.
#
#   ./failed_tasks.sh              # tally + failures of the last run
#   ./failed_tasks.sh <run_name>   # a named run (see: nextflow log)
#   ./failed_tasks.sh <run_name> all   # add the full per-task listing
#
# Reads the driver's own execution history, so no container is involved.
set -euo pipefail

run=${1:-last}
what=${2:-}

echo "== $run: tasks per process"
nextflow log "$run" -f process,status 2>/dev/null \
    | sort | uniq -c \
    | awk '{printf "  %-26s %-8s %s\n", $2, $3, $1}'

if [ "$what" = all ]; then
    echo
    echo "== all tasks"
    nextflow log "$run" -f name,status,exit,duration,workdir
fi

out=$(nextflow log "$run" -f name,status,exit,duration,workdir -F "status=='FAILED'" 2>/dev/null || true)
if [ -z "$out" ]; then
    echo
    echo "no failed tasks"
    exit 0
fi

echo
echo "== failures"
echo "$out"

# the first failure's stderr is usually the whole story
first=$(echo "$out" | head -1 | awk -F'\t' '{print $NF}')
if [ -s "$first/.command.err" ]; then
    echo
    echo "stderr of $(echo "$out" | head -1 | cut -f1):"
    tail -20 "$first/.command.err" | sed 's/^/  /'
fi
