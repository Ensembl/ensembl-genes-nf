#!/usr/bin/env bash
set -euo pipefail

bam=''; sample=''; mode=''; exclude_contigs=''; bed=''; bedgraph=''
while (($#)); do
    case "$1" in
        --bam) bam="$2"; shift 2 ;;
        --sample) sample="$2"; shift 2 ;;
        --mode) mode="$2"; shift 2 ;;
        --exclude-contigs) exclude_contigs="$2"; shift 2 ;;
        --bed) bed="$2"; shift 2 ;;
        --bedgraph) bedgraph="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done
[[ -n "$bam" && -n "$sample" && -n "$mode" && -n "$exclude_contigs" && -n "$bed" && -n "$bedgraph" ]] || { echo 'Missing CTSS argument' >&2; exit 2; }
case "$mode" in forward_read_5prime|reverse_read_5prime) ;; *) echo "Unsupported CTSS mode: $mode" >&2; exit 2 ;; esac
samtools view -F 4 "$bam" | awk -v sample="$sample" -v mode="$mode" -v excluded="$exclude_contigs" 'BEGIN {OFS="\t"; while ((getline line < excluded) > 0) blocked[line]=1} !blocked[$3] {chrom=$3; start=$4-1; flag=$2; read_length=length($10); reverse=int(flag/16)%2; if (mode == "reverse_read_5prime") reverse=1-reverse; if (reverse) {pos=start+read_length-1; strand="-"} else {pos=start; strand="+"}; key=chrom SUBSEP pos SUBSEP strand; count[key]++} END {for (key in count) {split(key, fields, SUBSEP); print fields[1], fields[2], fields[2]+1, sample, count[key], fields[3]}}' | sort -k1,1 -k2,2n > "$bed"
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$5}' "$bed" > "$bedgraph"
