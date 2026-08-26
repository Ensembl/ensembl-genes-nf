#!/usr/bin/env python3
"""Create deterministic transcript or genome partitions with workload estimates."""
import argparse
import csv
import os
import subprocess
import tempfile


def bed_records(path):
    records = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                records.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    return records


def gtf_records(path):
    records = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() not in {"gene", "transcript", "exon", "cds"}:
                continue
            marker = 'transcript_id "'
            if marker not in fields[8]:
                continue
            tid = fields[8].split(marker, 1)[1].split('"', 1)[0]
            start, end = int(fields[3]) - 1, int(fields[4])
            if tid not in records:
                records[tid] = [fields[0], start, end, tid]
            else:
                records[tid][1] = min(records[tid][1], start)
                records[tid][2] = max(records[tid][2], end)
    return list(records.values())


def fai_records(path):
    with open(path) as handle:
        return [(row[0], int(row[1])) for row in (line.rstrip("\n").split("\t") for line in handle) if len(row) >= 2]


def read_count(bam, regions):
    if not bam:
        return ""
    fd, region_file = tempfile.mkstemp(suffix=".bed")
    os.close(fd)
    try:
        with open(region_file, "w") as handle:
            for chrom, start, end, _ in regions:
                handle.write(f"{chrom}\t{start}\t{end}\n")
        result = subprocess.run(["samtools", "view", "-c", "-L", region_file, bam], check=True, capture_output=True, text=True)
        return int(result.stdout.strip())
    finally:
        os.unlink(region_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("transcriptome", "genome"), required=True)
    parser.add_argument("--bed12")
    parser.add_argument("--gtf")
    parser.add_argument("--fai")
    parser.add_argument("--bam")
    parser.add_argument("--partitions", type=int, default=1)
    parser.add_argument("--window-size", type=int, default=0)
    parser.add_argument("--padding", type=int, default=0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.partitions < 1 or args.padding < 0:
        parser.error("partitions must be positive and padding cannot be negative")
    if args.mode == "transcriptome":
        records = gtf_records(args.gtf) if args.gtf else bed_records(args.bed12)
        items = [(c, s, e, tid, 1) for c, s, e, tid in records]
        items.sort(key=lambda row: (row[0], row[1], row[3]))
    else:
        if not args.fai:
            parser.error("--fai is required for genome mode")
        window = args.window_size or 5000000
        items = [(c, s, min(s + window, length), f"{c}:{s}-{min(s + window, length)}", 0)
                 for c, length in fai_records(args.fai) for s in range(0, length, window)]
    count = min(args.partitions, len(items)) if items else 0
    buckets = [[] for _ in range(count)]
    for index, item in enumerate(items):
        buckets[index % count].append(item)
    assignments = {id(item): number for number, bucket in enumerate(buckets) for item in bucket}
    bucket_loads = [sum(row[4] for row in bucket) for bucket in buckets]
    bucket_reads = [read_count(args.bam, [(row[0], max(0, row[1] - args.padding), row[2] + args.padding, row[3]) for row in bucket]) if args.bam else "" for bucket in buckets]
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["partition_id", "mode", "contig", "start", "end", "padding", "annotation_load", "estimated_read_load"])
        for row in items:
            number = assignments[id(row)]
            writer.writerow([number, args.mode, row[0], max(0, row[1] - args.padding), row[2] + args.padding,
                             args.padding, bucket_loads[number], bucket_reads[number]])


if __name__ == "__main__":
    main()
