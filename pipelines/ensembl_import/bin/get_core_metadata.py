#!/usr/bin/env python3
"""Read the genebuild metadata required after importing a core database."""

from __future__ import annotations

import argparse
import json
import re

from _mysq_helper import mysql_connection

REQUIRED_METADATA = ("genebuild.version", "genebuild.last_geneset_update")
MONTH_PATTERN = re.compile(r"\d{4}-\d{2}")


def get_core_metadata(
    database: str,
    host: str,
    port: int,
    user: str,
    password: str,
) -> dict[str, str]:
    """Read and validate the genebuild version and update month."""
    with mysql_connection(database, host, port, user, password) as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT meta_key, meta_value
              FROM meta
             WHERE meta_key IN (%s, %s)
            """,
            REQUIRED_METADATA,
        )
        metadata = {key: value for key, value in cursor.fetchall()}

    missing = [key for key in REQUIRED_METADATA if not metadata.get(key)]
    if missing:
        raise ValueError(f"Missing required core metadata: {', '.join(missing)}")

    last_update = str(metadata["genebuild.last_geneset_update"])
    if not MONTH_PATTERN.fullmatch(last_update) or not 1 <= int(last_update[-2:]) <= 12:
        raise ValueError(
            "genebuild.last_geneset_update must be a valid month in YYYY-MM format; " f"found {last_update!r}"
        )

    return {
        "genebuild_version": str(metadata["genebuild.version"]),
        "last_geneset_update": last_update,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True)
    parser.add_argument("--host", required=True)
    parser.add_argument("--port", required=True, type=int)
    parser.add_argument("--user", required=True)
    parser.add_argument("--password", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(json.dumps(get_core_metadata(args.database, args.host, args.port, args.user, args.password)))


if __name__ == "__main__":
    main()
