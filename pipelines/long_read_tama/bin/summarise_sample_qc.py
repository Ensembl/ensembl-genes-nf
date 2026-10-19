#!/usr/bin/env python3
"""Join one-row-per-sample QC tables while preserving all source columns."""
import csv
import sys

sample, output, *inputs = sys.argv[1:]
merged = {"sample_id": sample}
for filename in inputs:
    with open(filename, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        continue
    row = rows[0]
    for key, value in row.items():
        if key and key != "sample_id":
            merged[key] = value

with open(output, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(merged), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerow(merged)
