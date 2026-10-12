#!/usr/bin/env python3
"""Calculate genome-wide annotation metrics from A's normalized TSV."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import heapq

from qc_utils import parse_biotypes, read_normalized, resolve_selected_features, write_metrics


def merge_intervals(features):
    by_seqid = defaultdict(list)
    for feature in features:
        by_seqid[feature.seqid].append((feature.start, feature.end))
    merged = {}
    covered = 0
    for seqid, intervals in by_seqid.items():
        result = []
        for start, end in sorted(intervals):
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], end))
            else:
                result.append((start, end))
        merged[seqid] = result
        covered += sum(end - start for start, end in result)
    return merged, covered


def overlap_counts(features):
    same = antisense = 0
    for seqid in sorted({x.seqid for x in features}):
        active = []
        by_end = []
        for feature in sorted((x for x in features if x.seqid == seqid), key=lambda x: (x.start, x.end)):
            while by_end and by_end[0][0] <= feature.start:
                _, index = heapq.heappop(by_end)
                active[index] = None
            for other in active:
                if other is None or other.end <= feature.start:
                    continue
                if other.strand == feature.strand:
                    same += 1
                elif {other.strand, feature.strand} == {"+", "-"}:
                    antisense += 1
            active.append(feature)
            heapq.heappush(by_end, (feature.end, len(active) - 1))
    return same, antisense


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("normalized")
    parser.add_argument("--feature-type", default="gene")
    parser.add_argument("--genome-size", type=int)
    parser.add_argument("--fasta-metrics")
    parser.add_argument("--biotype", action="append")
    parser.add_argument("--exclude-biotype", action="append")
    parser.add_argument("--canonical-only", action="store_true")
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--feature-overlaps-output")
    args = parser.parse_args()

    features = resolve_selected_features(
        read_normalized(args.normalized), parse_biotypes(args.biotype),
        parse_biotypes(args.exclude_biotype), args.canonical_only,
    )
    selected = [x for x in features if x.feature_type == args.feature_type]
    if not selected:
        parser.error(f"No {args.feature_type} features found")
    merged, covered = merge_intervals(selected)
    same, antisense = overlap_counts(selected)
    genome_size = args.genome_size
    if genome_size is None and args.fasta_metrics:
        with open(args.fasta_metrics, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row["metrics_name"] == "Genome size (bp)":
                    genome_size = int(row["metrics_value"])
                    break
    if genome_size is None:
        parser.error("--genome-size is required here; provide it from fasta_metrics.py")
    rows = [
        ("Genome size (bp)", genome_size),
        (f"{args.feature_type} count", len(selected)),
        (f"{args.feature_type} density ({args.feature_type}s/Mb)", f"{len(selected) / genome_size * 1e6:.4f}"),
        (f"% genome covered by {args.feature_type}", f"{covered / genome_size * 100:.4f}"),
        (f"Overlapping {args.feature_type} pairs (same strand)", same),
        (f"Antisense {args.feature_type} pairs (opposite strand)", antisense),
    ]
    if args.feature_overlaps_output:
        with open(args.feature_overlaps_output, "w", encoding="utf-8") as handle:
            handle.write("#chrom\toverlap_start\toverlap_end\tfeature1_id\tfeature1_start\tfeature1_end\tfeature1_strand\tfeature2_id\tfeature2_start\tfeature2_end\tfeature2_strand\trelation\toverlap_bp\n")
            for seqid in sorted({feature.seqid for feature in selected}):
                active = []
                for feature in sorted(
                    (x for x in selected if x.seqid == seqid),
                    key=lambda x: (x.start, x.end, x.feature_id),
                ):
                    active = [other for other in active if other.end > feature.start]
                    for other in active:
                        overlap_start = max(other.start, feature.start)
                        overlap_end = min(other.end, feature.end)
                        if overlap_end <= overlap_start:
                            continue
                        if other.strand == feature.strand:
                            relation = "same_strand"
                        elif {other.strand, feature.strand} == {"+", "-"}:
                            relation = "antisense"
                        else:
                            relation = "other_strand"
                        handle.write(
                            f"{seqid}\t{overlap_start}\t{overlap_end}\t"
                            f"{other.feature_id}\t{other.start}\t{other.end}\t{other.strand}\t"
                            f"{feature.feature_id}\t{feature.start}\t{feature.end}\t{feature.strand}\t"
                            f"{relation}\t{overlap_end - overlap_start}\n"
                        )
                    active.append(feature)
    write_metrics(rows, args.output_csv)


if __name__ == "__main__":
    main()

