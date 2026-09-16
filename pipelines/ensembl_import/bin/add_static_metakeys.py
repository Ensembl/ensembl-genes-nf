"""Add predefined metakeys to the core database from a static file."""

import argparse
import json

from _mysq_helper import mysql_connection


def parse_json(json_file: str) -> dict:
    """Parse a JSON file and return the data as a dictionary."""
    with open(json_file, "r", encoding="utf-8") as f:
        return json.load(f)


def add_metakeys(
    json_file: str, source: str, db_name: str, db_user: str, db_password: str, db_host: str, db_port: int
):
    """Add predefined metakeys to the core database from a static file."""
    data = parse_json(json_file)

    if source not in data:
        raise ValueError(f"'{source}' not found in {json_file}")

    metakeys = data[source]

    with mysql_connection(db_name, db_host, db_port, db_user, db_password) as conn:
        with conn.cursor() as cursor:
            for meta_key, meta_value in metakeys.items():
                cursor.execute(
                    """
                    INSERT IGNORE INTO meta (species_id, meta_key, meta_value)
                    VALUES (%s, %s, %s)
                    """,
                    (1, meta_key, meta_value),
                )
        conn.commit()

    print(f"Added {len(metakeys)} metakeys for '{source}'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json_file", type=str, help="Path to the JSON file")
    parser.add_argument("--source", type=str, help="Top-level JSON key to load, e.g. refseq")
    parser.add_argument("--db_name", type=str, help="Database name")
    parser.add_argument("--db_user", type=str, help="Database user")
    parser.add_argument("--db_password", type=str, help="Database password")
    parser.add_argument("--db_host", type=str, help="Database host")
    parser.add_argument("--db_port", type=int, help="Database port")

    args = parser.parse_args()

    add_metakeys(
        args.json_file,
        args.source,
        args.db_name,
        args.db_user,
        args.db_password,
        args.db_host,
        args.db_port,
    )
