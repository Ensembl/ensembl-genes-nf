#!/usr/bin/env python3
"""Build a sparse, wide junction matrix with two support columns per run."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_inputs(values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        run_id, separator, raw_path = value.partition("=")
        if not separator or not run_id or not raw_path:
            raise ValueError(f"invalid --input value: {value!r}")
        if run_id in result:
            raise ValueError(f"duplicate run ID: {run_id}")
        result[run_id] = Path(raw_path)
    if not result:
        raise ValueError("at least one input junction BED is required")
    return result


def build_matrix(inputs: dict[str, Path], output: Path, metadata_output: Path) -> None:
    records: dict[tuple[str, int, int, str], dict[str, object]] = {}
    for run_id, path in inputs.items():
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}:{line_number}: expected nine BED fields")
            chrom, start, end, name, _score, strand = fields[:6]
            unique_reads, multi_reads = fields[6:8]
            motif = fields[8]
            key = (chrom, int(start), int(end), strand)
            record = records.setdefault(key, {"name": name, "motifs": set(), "runs": defaultdict(lambda: [0, 0])})
            record["motifs"].add(motif)
            record["runs"][run_id][0] += int(unique_reads)
            record["runs"][run_id][1] += int(multi_reads)

    run_ids = list(inputs)
    fields = ["chrom", "chromStart", "chromEnd", "name", "strand", "splice_motif"]
    fields.extend(f"{run_id}.unique_reads" for run_id in run_ids)
    fields.extend(f"{run_id}.multi_reads" for run_id in run_ids)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for key in sorted(records, key=lambda item: (item[0], item[1], item[2], item[3])):
            chrom, start, end, strand = key
            record = records[key]
            motifs = {motif for motif in record["motifs"] if motif != "."}
            if len(motifs) > 1:
                raise ValueError(f"conflicting splice motifs for {chrom}:{start}-{end}:{strand}")
            motif = next(iter(motifs), ".")
            row = [chrom, start, end, record["name"], strand, motif]
            row.extend(record["runs"][run_id][0] for run_id in run_ids)
            row.extend(record["runs"][run_id][1] for run_id in run_ids)
            writer.writerow(row)

    metadata_output.write_text(json.dumps({
        "run_ids": run_ids,
        "columns": fields,
        "row_count": len(records),
    }, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    args = parser.parse_args()
    build_matrix(parse_inputs(args.input), args.output, args.metadata_output)


if __name__ == "__main__":
    main()
