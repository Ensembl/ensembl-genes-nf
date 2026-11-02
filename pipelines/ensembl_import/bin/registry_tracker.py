#!/usr/bin/env python3
"""Track RefSeq import pipelines in the assembly metadata registry.

Operations are ``pending``, ``start`` and ``complete``.  ``start`` prints the
new ``import_status_id``; pass that id to ``complete`` so a later attempt can
never update the wrong history row.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from dataclasses import dataclass
from typing import Any

import pymysql
from _mysq_helper import mysql_connection

LOGGER = logging.getLogger(__name__)

REGISTRY_SERVER = "mysql-ens-genebuild-prod-1"
REGISTRY_PORT = 4527
REGISTRY_USER = "ensadmin"
REGISTRY_PASSWORD = "ensembl"
REGISTRY_DB = "refseq_assembly_test"
GENEBUILD_VERSION = "IMP01"
GCF_ACCESSION_PATTERN = re.compile(r"GCF_\d+\.\d+")


class RegistryTrackerError(RuntimeError):
    """Raised when a registry operation cannot be completed safely."""


class ExistingImportError(RegistryTrackerError):
    """Raised when an assembly already has an import-status row."""


@dataclass(frozen=True)
class RegistryConfig:
    host: str = REGISTRY_SERVER
    port: int = REGISTRY_PORT
    user: str = REGISTRY_USER
    password: str = REGISTRY_PASSWORD
    database: str = REGISTRY_DB


PENDING_QUERY = """
    SELECT ra.*, s.*
      FROM refseq_assembly AS ra
      LEFT JOIN import_status AS s
        ON ra.refseq_assembly_id = s.refseq_assembly_id
     WHERE ra.is_current = 1
       AND s.import_status IS NULL
"""


def _assembly_lookup(cursor: pymysql.cursors.Cursor, refseq_accession: str) -> dict[str, Any]:
    cursor.execute(
        """
        SELECT ra.refseq_assembly_id, ra.refseq_accession, ra.is_current,
               s.import_status_id, s.import_status, s.date_started,
               s.date_status_update, s.genebuilder
          FROM refseq_assembly AS ra
          LEFT JOIN import_status AS s
            ON ra.refseq_assembly_id = s.refseq_assembly_id
         WHERE ra.refseq_accession = %s
           AND ra.is_current = 1
        """,
        (refseq_accession,),
    )
    rows = cursor.fetchall()
    if not rows:
        raise RegistryTrackerError(f"Current assembly not found: {refseq_accession}")
    if len(rows) > 1:
        raise RegistryTrackerError(f"Assembly {refseq_accession} has multiple import-status rows")
    row = rows[0]
    if isinstance(row, dict):
        return row
    return dict(zip((column[0] for column in cursor.description), row))


def list_pending(config: RegistryConfig) -> list[dict[str, Any]]:
    """Return current RefSeq assemblies without an import-status record."""
    with mysql_connection(
        config.database, config.host, config.port, config.user, config.password
    ) as conn, conn.cursor(pymysql.cursors.DictCursor) as cursor:
        cursor.execute(PENDING_QUERY)
        return list(cursor.fetchall())


def _acquire_assembly_lock(cursor: pymysql.cursors.Cursor, assembly_id: int) -> str:
    """Serialize start operations for one assembly on this MySQL server."""
    lock_name = f"ensembl_import:{assembly_id}"
    cursor.execute("SELECT GET_LOCK(%s, 30)", (lock_name,))
    if cursor.fetchone()[0] != 1:
        raise RegistryTrackerError(f"Could not acquire registry lock for assembly {assembly_id}")
    return lock_name


def start_import(config: RegistryConfig, refseq_accession: str, genebuilder: str) -> int:
    """Create an in-progress row and return its id."""
    if not GCF_ACCESSION_PATTERN.fullmatch(refseq_accession):
        raise RegistryTrackerError(
            f"Invalid refseq_accession: {refseq_accession}; expected a GCF accession such as GCF_000001405.40"
        )
    if len(genebuilder) > 20:
        raise RegistryTrackerError("genebuilder must be at most 20 characters")

    with mysql_connection(
        config.database, config.host, config.port, config.user, config.password
    ) as conn, conn.cursor(pymysql.cursors.DictCursor) as cursor:
        assembly = _assembly_lookup(cursor, refseq_accession)
        refseq_accession = assembly["refseq_accession"]
        lock_name = _acquire_assembly_lock(cursor, assembly["refseq_assembly_id"])
        try:
            # Re-check after locking: another launcher may have inserted a
            # row between the first lookup and lock acquisition.
            cursor.execute(
                "SELECT import_status_id, import_status, genebuilder "
                "FROM import_status WHERE refseq_assembly_id = %s",
                (assembly["refseq_assembly_id"],),
            )
            existing = cursor.fetchall()
            if existing:
                details = ", ".join(
                    f"id={row['import_status_id']} status={row['import_status']} builder={row['genebuilder']}"
                    for row in existing
                )
                raise ExistingImportError(
                    f"Import already recorded for {refseq_accession}; refusing to start another ({details})"
                )

            cursor.execute(
                """
                INSERT INTO import_status
                    (refseq_assembly_id, import_status, genebuilder,
                     annotation_source, date_started, date_status_update,
                     genebuild_version)
                VALUES (%s, 'in_progress', %s, 'import_refseq',
                        CURRENT_DATE, CURRENT_DATE, %s)
                """,
                (assembly["refseq_assembly_id"], genebuilder, GENEBUILD_VERSION),
            )
            status_id = cursor.lastrowid
            conn.commit()
            LOGGER.info("Started import %s with import_status_id=%s", refseq_accession, status_id)
            return int(status_id)
        finally:
            cursor.execute("SELECT RELEASE_LOCK(%s)", (lock_name,))


def complete_import(config: RegistryConfig, import_status_id: int) -> None:
    """Mark one in-progress import as completed."""
    with mysql_connection(
        config.database, config.host, config.port, config.user, config.password
    ) as conn, conn.cursor() as cursor:
        cursor.execute(
            """
            UPDATE import_status
               SET import_status = 'completed',
                   date_status_update = CURRENT_DATE,
                   last_genebuild_update = CURRENT_DATE,
                   last_attempt = 1
             WHERE import_status_id = %s
               AND import_status = 'in_progress'
            """,
            (import_status_id,),
        )
        if cursor.rowcount != 1:
            conn.rollback()
            raise RegistryTrackerError(f"No in-progress import found for import_status_id={import_status_id}")
        conn.commit()
    LOGGER.info("Completed import_status_id=%s", import_status_id)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=REGISTRY_SERVER)
    parser.add_argument("--port", default=REGISTRY_PORT, type=int)
    parser.add_argument("--user", default=REGISTRY_USER)
    parser.add_argument("--password", default=REGISTRY_PASSWORD)
    parser.add_argument("--database", default=REGISTRY_DB)
    subparsers = parser.add_subparsers(dest="operation", required=True)
    subparsers.add_parser("pending", help="list current assemblies with no import status")
    start = subparsers.add_parser("start", help="record an import start")
    start.add_argument("--accession", required=True, dest="refseq_accession", help="current GCF accession")
    start.add_argument("--genebuilder", required=True, help="registry genebuilder")
    complete = subparsers.add_parser("complete", help="record a successful import")
    complete.add_argument("--import-status-id", required=True, type=int)
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    config = RegistryConfig(args.host, args.port, args.user, args.password, args.database)
    if args.operation == "pending":
        print(json.dumps(list_pending(config), default=str, indent=2))
    elif args.operation == "start":
        print(start_import(config, args.refseq_accession, args.genebuilder))
    else:
        complete_import(config, args.import_status_id)


if __name__ == "__main__":
    main()
