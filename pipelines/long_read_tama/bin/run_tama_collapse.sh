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
tama_args=("$@")

status_file="tama_status.tsv"
stderr_file="tama_collapse.stderr"

run_collapse() {
    local input_bam="$1"
    local stderr_out="$2"
    tama_collapse.py -s "$input_bam" -b BAM -f "$reference" -p "$prefix" "${tama_args[@]}" \
        > tama_collapse.stdout 2> "$stderr_out" || return $?
    return 0
}

if run_collapse "$bam" "$stderr_file"; then
    rc=0
else
    rc=$?
fi

# TAMA 1.0.3 can leave an empty internal locus after multi-map resolution and
# then dereference trans_list[0]. Retry this specific upstream bug with only
# primary mapped alignments. Do not apply the workaround to unrelated errors.
if [[ "$rc" -ne 0 ]] && grep -Fq 'IndexError: list index out of range' "$stderr_file"; then
    filtered_bam="${prefix}.primary.bam"
    retry_stderr="${prefix}.retry.stderr"
    if samtools view -bh -F 2308 -o "$filtered_bam" "$bam" \
            && samtools index "$filtered_bam" \
            && run_collapse "$filtered_bam" "$retry_stderr"; then
        printf 'TAMA workaround: retried after filtering to primary mapped alignments\n' > "$stderr_file"
        cat "$retry_stderr" >> "$stderr_file"
        rc=0
    else
        retry_rc=$?
        cat "$retry_stderr" >> "$stderr_file" 2>/dev/null || true
        rc="$retry_rc"
    fi
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
