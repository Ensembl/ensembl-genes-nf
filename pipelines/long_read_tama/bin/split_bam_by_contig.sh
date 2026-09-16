#!/usr/bin/env bash
set -euo pipefail
bam=$1; workload=$2; accession=$3; shard_dir=$4; manifest=$5
mkdir -p "$shard_dir"
printf 'contig\treference_bases\tmapped_reads\tresource_class\tpath\n' > "$manifest"
tail -n +2 "$workload" | while IFS=$(printf '\t') read -r contig length mapped eligible resource_class; do
    safe=$(printf '%s' "$contig" | tr -c 'A-Za-z0-9._-' '_')
    shard="$shard_dir/$accession.$safe.bam"
    samtools view -b -o "$shard" "$bam" "$contig"
    samtools index "$shard"
    printf '%s\t%s\t%s\t%s\t%s\n' "$contig" "$length" "$mapped" "$resource_class" "$(basename "$shard")" >> "$manifest"
done
test "$(wc -l < "$manifest")" -gt 1 || { echo 'No contig shards were materialised' >&2; exit 1; }
