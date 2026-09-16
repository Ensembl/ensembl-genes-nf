#!/usr/bin/env bash
set -euo pipefail
bam=$1; bai=$2; output=$3; threshold=$4
test -r "$bam" && test -r "$bai" || { echo 'BAM and BAI must be readable' >&2; exit 1; }
samtools quickcheck -v "$bam"
samtools idxstats "$bam" | awk -v threshold="$threshold" 'BEGIN { OFS="\t"; print "contig", "reference_bases", "mapped_reads", "eligible", "resource_class" }
    $1 != "*" && $2 > 0 { resource = ($3 >= threshold * 5 ? "very_large" : ($3 >= threshold ? "large" : "small")); print $1, $2, $3, "true", resource }' > "$output"
test "$(wc -l < "$output")" -gt 1 || { echo 'BAM has no eligible reference contigs' >&2; exit 1; }
