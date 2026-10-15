#!/usr/bin/env python3
"""Scan FASTA once for genome size, masking, and exon-GC metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict

from qc_utils import parse_biotypes, read_normalized, resolve_selected_features, write_metrics


LOWERCASE_RUN = re.compile(rb"[a-z]+")


def intervals_by_seqid(features):
    result = defaultdict(list)
    for index, feature in enumerate(features):
        result[feature.seqid].append((feature.start, feature.end, index, feature))
    for values in result.values():
        values.sort()
    return result


def fasta_seqid(header):
    return header[1:].split(None, 1)[0].decode()


def scan(path, exons, annotations, masked_output=None):
    exon_intervals = intervals_by_seqid(exons)
    annotation_intervals = intervals_by_seqid(annotations)
    gc_counts = [0] * len(exons)
    base_counts = [0] * len(exons)
    exon_cursors = {}
    annotation_cursors = {}
    genome_size = total_masked = masked_annotated = 0
    seen = set()
    overlap_out = None
    if masked_output:
        overlap_out = gzip.open(masked_output, "wt", encoding="utf-8") if str(masked_output).endswith(".gz") else open(masked_output, "wt", encoding="utf-8")
        overlap_out.write("#chrom\toverlap_start\toverlap_end\tmask_start\tmask_end\tfeature_id\tfeature_start\tfeature_end\tfeature_strand\toverlap_bp\n")

    try:
        handle = gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")
        with handle:
            seqid = None
            position = 0
            mask_start = None
            mask_end = None

            def finish_mask():
                nonlocal total_masked, masked_annotated, mask_start, mask_end
                if seqid is None or mask_start is None or mask_end is None:
                    return
                mask_length = mask_end - mask_start
                total_masked += mask_length
                annotated = False
                intervals = annotation_intervals.get(seqid, ())
                cursor = annotation_cursors.get(seqid, 0)
                while cursor < len(intervals) and intervals[cursor][1] <= mask_start:
                    cursor += 1
                for interval_index in range(cursor, len(intervals)):
                    start, end, _index, feature = intervals[interval_index]
                    if start >= mask_end:
                        break
                    overlap_start = max(mask_start, start)
                    overlap_end = min(mask_end, end)
                    if overlap_end <= overlap_start:
                        continue
                    annotated = True
                    if overlap_out:
                        overlap_out.write(f"{seqid}\t{overlap_start}\t{overlap_end}\t{mask_start}\t{mask_end}\t{feature.feature_id}\t{feature.start}\t{feature.end}\t{feature.strand}\t{overlap_end - overlap_start}\n")
                while cursor < len(intervals) and intervals[cursor][1] <= mask_end:
                    cursor += 1
                annotation_cursors[seqid] = cursor
                if annotated:
                    masked_annotated += mask_length
                mask_start = mask_end = None

            for raw_line in handle:
                if raw_line.startswith(b">"):
                    finish_mask()
                    seqid = fasta_seqid(raw_line)
                    seen.add(seqid)
                    position = 0
                    continue
                sequence = raw_line.strip()
                if not sequence or seqid is None:
                    continue
                line_start, line_end = position, position + len(sequence)
                genome_size += len(sequence)
                for match in LOWERCASE_RUN.finditer(sequence):
                    start, end = line_start + match.start(), line_start + match.end()
                    if mask_start is not None:
                        if start == mask_end:
                            mask_end = end
                            continue
                        finish_mask()
                    mask_start, mask_end = start, end
                if mask_start is not None and mask_end != line_end:
                    finish_mask()

                intervals = exon_intervals.get(seqid, ())
                cursor = exon_cursors.get(seqid, 0)
                while cursor < len(intervals) and intervals[cursor][1] <= line_start:
                    cursor += 1
                for interval_index in range(cursor, len(intervals)):
                    start, end, index, _feature = intervals[interval_index]
                    if start >= line_end:
                        break
                    overlap_start, overlap_end = max(line_start, start), min(line_end, end)
                    if overlap_end <= overlap_start:
                        continue
                    chunk = sequence[overlap_start - line_start:overlap_end - line_start].upper()
                    gc_counts[index] += chunk.count(b"G") + chunk.count(b"C")
                    base_counts[index] += sum(chunk.count(base) for base in b"ACGT")
                while cursor < len(intervals) and intervals[cursor][1] <= line_end:
                    cursor += 1
                exon_cursors[seqid] = cursor
                position = line_end
            finish_mask()
    finally:
        if overlap_out:
            overlap_out.close()

    gc_values = [(gc / bases) * 100 for gc, bases in zip(gc_counts, base_counts) if bases]
    return genome_size, total_masked, masked_annotated, gc_values, len(gc_values), sorted(set(exon_intervals) - seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("normalized")
    parser.add_argument("fasta")
    parser.add_argument("--feature-type", default="gene")
    parser.add_argument("--biotype", action="append")
    parser.add_argument("--exclude-biotype", action="append")
    parser.add_argument("--canonical-only", action="store_true")
    parser.add_argument("--gc-bin-size", type=int, default=10)
    parser.add_argument("--masked-overlaps-output")
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()
    if not 1 <= args.gc_bin_size <= 100:
        parser.error("--gc-bin-size must be between 1 and 100")

    features = resolve_selected_features(read_normalized(args.normalized), parse_biotypes(args.biotype), parse_biotypes(args.exclude_biotype), args.canonical_only)
    exons = [feature for feature in features if feature.feature_type.lower() == "exon"]
    annotations = [feature for feature in features if feature.feature_type == args.feature_type]
    genome_size, masked, masked_annotated, gc_values, gc_count, missing = scan(args.fasta, exons, annotations, args.masked_overlaps_output)
    rows = [("Genome size (bp)", genome_size), ("Total soft-masked bases", f"{masked} ({masked / genome_size * 100:.2f}%)" if genome_size else masked), ("% masked regions annotated", "N/A (no soft-masking detected)" if not masked else f"{masked_annotated / masked * 100:.4f}"), ("Exon GC count", gc_count), ("Exons skipped for GC", len(exons) - gc_count)]
    if gc_values:
        ordered_gc = sorted(gc_values)
        midpoint = len(ordered_gc) // 2
        gc_median = ordered_gc[midpoint] if len(ordered_gc) % 2 else (ordered_gc[midpoint - 1] + ordered_gc[midpoint]) / 2
        rows.extend([("Exon GC mean (%)", f"{sum(gc_values) / len(gc_values):.2f}"), ("Exon GC median (%)", f"{gc_median:.2f}"), ("Exon GC min (%)", f"{min(gc_values):.2f}"), ("Exon GC max (%)", f"{max(gc_values):.2f}")])
    else:
        rows.extend((f"Exon GC {label} (%)", "N/A") for label in ("mean", "median", "min", "max"))
    if missing:
        rows.extend([("Annotation seqids missing from FASTA", len(missing)), ("Missing annotation seqids", ", ".join(missing[:10]))])
    for index in range((100 + args.gc_bin_size - 1) // args.gc_bin_size):
        start, end = index * args.gc_bin_size, min(100, (index + 1) * args.gc_bin_size)
        count = sum(start <= value < end if end < 100 else start <= value <= end for value in gc_values)
        label = f"{start}-<{end}%" if end < 100 else f"{start}-100%"
        rows.append((f"Exon GC {label}", f"{count} ({count / len(gc_values) * 100:.2f}%)" if gc_values else "0 (N/A)"))
    write_metrics(rows, args.output_csv)


if __name__ == "__main__":
    main()
