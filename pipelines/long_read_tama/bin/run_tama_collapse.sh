#!/usr/bin/env bash

set -euo pipefail

if (( $# < 6 )); then
    echo "usage: run_tama_collapse.sh ACCESSION SHARD BAM REFERENCE PREFIX ARGS..." >&2
    exit 2
fi

accession="$1"
shard="$2"
bam="$3"
reference="$4"
prefix="$5"
shift 5

status_file="tama_status.tsv"
stderr_file="tama_collapse.stderr"

if tama_collapse.py -s "$bam" -b BAM -f "$reference" -p "$prefix" "$@" \
        > tama_collapse.stdout 2> "$stderr_file"; then
    rc=0
else
    rc=$?
fi

if [[ "$rc" -eq 0 && -s "${prefix}_collapsed.bed" ]]; then
    status="SUCCESS"
else
    status="TAMA_FAILED"
    rm -f "${prefix}_collapsed.bed" "${prefix}_read.txt" "${prefix}_trans_read.bed"
fi

printf 'accession\tshard\tstatus\texit_status\n%s\t%s\t%s\t%s\n' \
    "$accession" "$shard" "$status" "$rc" > "$status_file"

# A shard-level TAMA failure is data recorded for the cohort, not a workflow
# failure. The status and stderr outputs are the authoritative diagnostics.
exit 0
