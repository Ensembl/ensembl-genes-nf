#!/usr/bin/env python
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Utilities for statistics pipelines."""

# from __future__ import annotations

import argparse
from typing import List
import pymysql


def get_meta_value(
    *,
    dbname: str,
    meta_key: str,
    species_id: int = 1,
    host: str,
    port: int,
    user: str,
    password: str,
) -> List[str]:
    """
    Retrieve meta_value(s) from an Ensembl core database.

    Args:
        dbname: Core database name
        meta_key: Meta key to query
        species_id: Species ID (default: 1)
        host: MySQL host
        port: MySQL port
        user: MySQL user
        password: MySQL password

    Returns:
        List of meta_value strings (possibly empty).
    """

    query = """
        SELECT meta_value
        FROM meta
        WHERE meta_key = %s
        AND species_id = %s
    """

    conn = pymysql.connect(
        host=host,
        port=port,
        user=user,
        password=password,
        database=dbname,
        cursorclass=pymysql.cursors.Cursor,
    )

    try:
        with conn.cursor() as cursor:
            cursor.execute(query, (meta_key, species_id))
            rows = cursor.fetchall()
            return [row[0] for row in rows]
    finally:
        conn.close()


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fetch meta_value from an Ensembl core database")
    parser.add_argument("--db", required=True, help="Core database name")
    parser.add_argument("--key", required=True, help="Meta key to fetch")
    parser.add_argument("--species-id", type=int, default=1)
    parser.add_argument("--host", required=True)
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--user", required=True)
    parser.add_argument("--password", default="")
    return parser


def main() -> None:
    """Main entry-point."""
    parser = _build_arg_parser()
    args = parser.parse_args()

    values = get_meta_value(
        dbname=args.db,
        meta_key=args.key,
        species_id=args.species_id,
        host=args.host,
        port=args.port,
        user=args.user,
        password=args.password,
    )

    for value in values:
        print(value)


if __name__ == "__main__":
    main()
