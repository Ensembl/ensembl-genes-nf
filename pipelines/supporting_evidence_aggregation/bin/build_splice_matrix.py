#!/usr/bin/env python3
"""Build a sparse, wide junction matrix with two support columns per run."""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
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
    run_ids = list(inputs)
    run_indices = {run_id: index for index, run_id in enumerate(run_ids)}
    fields = ["chrom", "chromStart", "chromEnd", "name", "strand", "splice_motif"]
    fields.extend(f"{run_id}.unique_reads" for run_id in run_ids)
    fields.extend(f"{run_id}.multi_reads" for run_id in run_ids)

    # The union of junctions across hundreds of runs can be much larger than
    # the task's memory allocation. Keep the sparse support matrix on disk and
    # stream the sorted join into the wide output table.
    database = output.with_name(f".{output.name}.sqlite")
    connection = sqlite3.connect(database)
    try:
        connection.executescript(
            """
            PRAGMA journal_mode=OFF;
            PRAGMA synchronous=OFF;
            PRAGMA temp_store=FILE;
            PRAGMA cache_size=-262144;
            CREATE TABLE junctions (
                chrom TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,
                strand TEXT NOT NULL,
                name TEXT NOT NULL,
                motif TEXT NOT NULL,
                PRIMARY KEY (chrom, start, end, strand)
            );
            CREATE TABLE support (
                chrom TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,
                strand TEXT NOT NULL,
                run_index INTEGER NOT NULL,
                unique_reads INTEGER NOT NULL,
                multi_reads INTEGER NOT NULL,
                PRIMARY KEY (chrom, start, end, strand, run_index)
            );
            """
        )
        junction_statement = """
            INSERT INTO junctions(chrom, start, end, strand, name, motif)
            VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT(chrom, start, end, strand) DO UPDATE SET
                motif = CASE
                    WHEN junctions.motif = '.' THEN excluded.motif
                    WHEN excluded.motif = '.' OR junctions.motif = excluded.motif THEN junctions.motif
                    ELSE '__CONFLICT__'
                END
        """
        support_statement = """
            INSERT INTO support(chrom, start, end, strand, run_index, unique_reads, multi_reads)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(chrom, start, end, strand, run_index) DO UPDATE SET
                unique_reads = support.unique_reads + excluded.unique_reads,
                multi_reads = support.multi_reads + excluded.multi_reads
        """
        for run_id, path in inputs.items():
            junction_rows = []
            support_rows = []
            with path.open(encoding="utf-8") as handle:
                for line_number, line in enumerate(handle, 1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    values = line.rstrip("\n").split("\t")
                    if len(values) < 9:
                        raise ValueError(f"{path}:{line_number}: expected nine BED fields")
                    chrom, raw_start, raw_end, name, _score, strand = values[:6]
                    unique_reads, multi_reads = values[6:8]
                    motif = values[8]
                    start = int(raw_start)
                    end = int(raw_end)
                    junction_rows.append((chrom, start, end, strand, name, motif))
                    support_rows.append((chrom, start, end, strand, run_indices[run_id], int(unique_reads), int(multi_reads)))
                    if len(junction_rows) >= 10000:
                        connection.executemany(junction_statement, junction_rows)
                        connection.executemany(support_statement, support_rows)
                        junction_rows.clear()
                        support_rows.clear()
            if junction_rows:
                connection.executemany(junction_statement, junction_rows)
                connection.executemany(support_statement, support_rows)
        connection.commit()

        conflict = connection.execute(
            "SELECT chrom, start, end, strand FROM junctions WHERE motif = '__CONFLICT__' LIMIT 1"
        ).fetchone()
        if conflict:
            chrom, start, end, strand = conflict
            raise ValueError(f"conflicting splice motifs for {chrom}:{start}-{end}:{strand}")

        row_count = 0
        current_key = None
        current_name = None
        current_motif = None
        unique_values = [0] * len(run_ids)
        multi_values = [0] * len(run_ids)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            cursor = connection.execute(
                """
                SELECT j.chrom, j.start, j.end, j.name, j.strand, j.motif,
                       s.run_index, s.unique_reads, s.multi_reads
                FROM junctions AS j
                JOIN support AS s
                  ON s.chrom = j.chrom AND s.start = j.start
                 AND s.end = j.end AND s.strand = j.strand
                ORDER BY j.chrom, j.start, j.end, j.strand, s.run_index
                """
            )
            for chrom, start, end, name, strand, motif, run_index, unique_reads, multi_reads in cursor:
                key = (chrom, start, end, strand)
                if current_key is not None and key != current_key:
                    writer.writerow([current_key[0], current_key[1], current_key[2], current_name, current_key[3], current_motif, *unique_values, *multi_values])
                    row_count += 1
                    unique_values = [0] * len(run_ids)
                    multi_values = [0] * len(run_ids)
                if key != current_key:
                    current_key = key
                    current_name = name
                    current_motif = motif
                unique_values[run_index] = unique_reads
                multi_values[run_index] = multi_reads
            if current_key is not None:
                writer.writerow([current_key[0], current_key[1], current_key[2], current_name, current_key[3], current_motif, *unique_values, *multi_values])
                row_count += 1
    finally:
        connection.close()
        database.unlink(missing_ok=True)

    metadata_output.write_text(json.dumps({
        "run_ids": run_ids,
        "columns": fields,
        "row_count": row_count,
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
