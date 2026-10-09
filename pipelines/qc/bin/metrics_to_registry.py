#!/usr/bin/env python3
"""Load the final AGAT metrics CSV into registry new_metrics."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger("metrics_to_registry")


def read_metrics(csv_path: str | Path) -> list[tuple[str, str]]:
    """Read the two-column metric-name/metric-value AGAT output."""

    metrics: list[tuple[str, str]] = []
    with Path(csv_path).open(newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle), start=1):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) < 2:
                raise ValueError(
                    f"{csv_path}:{line_number} has fewer than two CSV columns"
                )

            name = row[0].strip()
            value = row[1].strip()
            if not name:
                LOGGER.warning("Skipping unnamed metric on line %s", line_number)
                continue
            if line_number == 1 and name.lower() in {
                "metric",
                "metric_name",
                "name",
            }:
                continue
            metrics.append((name, value))

    if not metrics:
        raise ValueError(f"No metrics found in {csv_path}")
    return metrics


def _connect(args: argparse.Namespace) -> Any:
    """Create the registry connection lazily so parsing can be tested offline."""

    import pymysql

    return pymysql.connect(
        host=args.registry_host,
        port=args.registry_port,
        user=args.registry_user,
        password=args.registry_password,
        database=args.registry_db,
        cursorclass=pymysql.cursors.DictCursor,
    )


def fetch_registry_ids(
    connection: Any,
    assembly: str,
) -> tuple[int, int]:
    """Find the assembly and latest live genebuild status for an accession."""

    query = """
        SELECT a.assembly_id, gs.genebuild_status_id
        FROM assembly AS a
        JOIN genebuild_status AS gs ON gs.assembly_id = a.assembly_id
        WHERE gs.gca_accession = %s
          AND gs.gb_status = 'live'
        ORDER BY gs.genebuild_status_id DESC
        LIMIT 1
    """

    with connection.cursor() as cursor:
        cursor.execute(query, [assembly])
        row = cursor.fetchone()
    if not row:
        raise LookupError(f"No live registry genebuild found for {assembly}")
    return int(row["assembly_id"]), int(row["genebuild_status_id"])


def write_metrics(
    connection: Any,
    assembly_id: int,
    genebuild_status_id: int,
    metrics: list[tuple[str, str]],
) -> None:
    """Replace the supplied metric names atomically in new_metrics."""

    metric_names = [name for name, _value in metrics]
    placeholders = ",".join(["%s"] * len(metric_names))
    delete_query = f"""
        DELETE FROM new_metrics
        WHERE genebuild_status_id = %s
          AND metrics_name IN ({placeholders})
    """
    insert_query = """
        INSERT INTO new_metrics
            (genebuild_status_id, assembly_id, metrics_name, metrics_value)
        VALUES (%s, %s, %s, %s)
    """

    with connection.cursor() as cursor:
        cursor.execute(delete_query, [genebuild_status_id, *metric_names])
        cursor.executemany(
            insert_query,
            [
                (genebuild_status_id, assembly_id, metric_name, metric_value)
                for metric_name, metric_value in metrics
            ],
        )


def main(args: argparse.Namespace) -> None:
    """Load one AGAT metrics file into the registry."""

    metrics = read_metrics(args.metrics_csv)
    connection = _connect(args)
    try:
        assembly_id, genebuild_status_id = fetch_registry_ids(
            connection,
            args.assembly,
        )
        write_metrics(connection, assembly_id, genebuild_status_id, metrics)
        connection.commit()
        LOGGER.info(
            "Wrote %s metrics for %s (assembly_id=%s, genebuild_status_id=%s)",
            len(metrics),
            args.assembly,
            assembly_id,
            genebuild_status_id,
        )
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser used by the Nextflow module."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics_csv", help="Final AGAT parsed metrics CSV")
    parser.add_argument("--assembly", required=True, help="GCA assembly accession")
    parser.add_argument("--registry-host", required=True)
    parser.add_argument("--registry-port", type=int, required=True)
    parser.add_argument("--registry-user", required=True)
    parser.add_argument("--registry-password", required=True)
    parser.add_argument("--registry-db", required=True)
    return parser


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main(build_parser().parse_args())
