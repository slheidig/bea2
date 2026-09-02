#!/usr/bin/env bash
# Which tasks failed, and where to look. Replaces -with-trace, which cannot be
# used here: it injects a ps-based collector into each container and none of our
# images ship procps, so every containerised task would exit 1 immediately.
#
#   ./failed_tasks.sh              # the last run
#   ./failed_tasks.sh <run_name>   # a named run (see: nextflow log)
#   ./failed_tasks.sh <run_name> all   # every task, not only the failures
#
# Reads the driver's own execution history, so no container is involved.
set -euo pipefail

run=${1:-last}
what=${2:-failed}
fields='name,status,exit,duration,workdir'

if [ "$what" = all ]; then
    nextflow log "$run" -f "$fields"
    exit 0
fi

out=$(nextflow log "$run" -f "$fields" -F "status=='FAILED'")
if [ -z "$out" ]; then
    echo "no failed tasks in run '$run'"
    exit 0
fi

echo "$out"
echo
echo "affected processes:"
echo "$out" | awk -F'\t' '{p=$1; sub(/ .*/,"",p); c[p]++} END{for (k in c) printf "  %-22s %d\n", k, c[k]}' | sort

# the first failure's stderr is usually the whole story
first=$(echo "$out" | head -1 | awk -F'\t' '{print $NF}')
if [ -s "$first/.command.err" ]; then
    echo
    echo "stderr of $(echo "$out" | head -1 | cut -f1):"
    tail -20 "$first/.command.err" | sed 's/^/  /'
fi
