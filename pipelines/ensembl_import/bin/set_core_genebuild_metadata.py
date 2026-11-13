#!/usr/bin/env python3
"""Write and return the genebuild metadata for a newly loaded core."""

from __future__ import annotations

import argparse
import json
from datetime import date

from _mysq_helper import mysql_connection

MONTH_FORMAT = "%Y-%m"


def set_core_genebuild_metadata(
    database: str,
    host: str,
    port: int,
    user: str,
    password: str,
    genebuild_version: str,
) -> dict[str, str]:
    """Write the metadata and return the exact values written."""
    last_geneset_update = date.today().strftime(MONTH_FORMAT)

    with mysql_connection(database, host, port, user, password) as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            INSERT INTO meta (species_id, meta_key, meta_value)
            VALUES
                (1, 'genebuild.version', %s),
                (1, 'genebuild.last_geneset_update', %s)
            ON DUPLICATE KEY UPDATE meta_value = VALUES(meta_value)
            """,
            (genebuild_version, last_geneset_update),
        )
        connection.commit()

    return {
        "genebuild_version": genebuild_version,
        "last_geneset_update": last_geneset_update,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True)
    parser.add_argument("--host", required=True)
    parser.add_argument("--port", required=True, type=int)
    parser.add_argument("--user", required=True)
    parser.add_argument("--password", default="")
    parser.add_argument("--genebuild-version", required=True)
    args = parser.parse_args()

    print(
        json.dumps(
            set_core_genebuild_metadata(
                args.database,
                args.host,
                args.port,
                args.user,
                args.password,
                args.genebuild_version,
            )
        )
    )


if __name__ == "__main__":
    main()
