#!/usr/bin/env python3
"""Summarise Nextflow trace rows by candidate-model backend."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


BACKENDS = ("tama", "stringtie2", "stringtie3", "tmerge", "isoquant", "flair", "bambu")


def backend_for(name: str) -> str:
    lowered = name.lower()
    for backend in BACKENDS:
        if backend in lowered:
            return backend
    return "shared"


def seconds(value: str) -> float:
    value = value.strip()
    if not value:
        return 0.0
    total = 0.0
    for number, unit in re.findall(r"([0-9]+(?:\.[0-9]+)?)\s*([dhms])", value.lower()):
        total += float(number) * {"d": 86400, "h": 3600, "m": 60, "s": 1}[unit]
    if total:
        return total
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*ms", value.lower())
    return float(match.group(1)) / 1000 if match else 0.0


def gigabytes(value: str) -> float:
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*([kmgt]?b)", value.lower())
    if not match:
        return 0.0
    number = float(match.group(1))
    scale = {"b": 1 / 1024**3, "kb": 1 / 1024**2, "mb": 1 / 1024, "gb": 1, "tb": 1024}
    return number * scale[match.group(2)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trace", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    groups: dict[str, dict[str, object]] = defaultdict(lambda: {
        "task_count": 0, "completed_count": 0, "failed_count": 0,
        "ignored_count": 0, "total_realtime_seconds": 0.0, "peak_rss_gb": 0.0,
    })
    with args.trace.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            group = groups[backend_for(row.get("name", ""))]
            group["task_count"] += 1
            status = row.get("status", "").upper()
            group[f"{status.lower()}_count"] = group.get(f"{status.lower()}_count", 0) + 1
            group["total_realtime_seconds"] += seconds(row.get("realtime", ""))
            group["peak_rss_gb"] = max(float(group["peak_rss_gb"]), gigabytes(row.get("peak_rss", "")))

    rows = [{"backend": backend, **values} for backend, values in sorted(groups.items())]
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    fields = ["backend", "task_count", "completed_count", "failed_count", "ignored_count", "total_realtime_seconds", "peak_rss_gb"]
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(json.dumps({"groups": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
