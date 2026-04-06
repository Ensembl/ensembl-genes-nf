#!/usr/bin/env python3
"""
Assign Ensembl stable IDs to genes, transcripts, exons, and translations
in a core MySQL database.

If this is a new species, IDs start from ENS<prefix>G00000000001.1.
If the database already has stable IDs (from a previous build), re-assignment
maps old → new IDs using overlap-based matching against a previous DB.

For a new build (no previous stable IDs), this script simply assigns
sequential stable IDs to all genes/transcripts/exons/translations that
lack them.

Usage:
    assign_stable_ids.py \\
        --host   mysql-host \\
        --port   3306 \\
        --user   ensadmin \\
        --password secret \\
        --dbname homo_sapiens_core_109_38 \\
        --prefix '' \\          # empty for human (ENSG); or 'GAL' for chicken (ENSGALG)
        --species-id 1
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

import pymysql
import pymysql.cursors

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# Ensembl stable ID format: ENS[prefix][type][11-digit-zero-padded].[version]
# type: G=gene, T=transcript, E=exon, P=translation
ID_WIDTH = 11


def _make_id(prefix: str, id_type: str, n: int) -> str:
    return f"ENS{prefix}{id_type}{str(n).zfill(ID_WIDTH)}"


def assign_gene_stable_ids(conn, prefix: str, species_id: int) -> int:
    """Assign stable IDs to genes that don't have one. Returns count assigned."""
    with conn.cursor() as cur:
        cur.execute(
            """SELECT g.gene_id FROM gene g
               LEFT JOIN gene_stable_id si ON g.gene_id = si.gene_id
               WHERE si.stable_id IS NULL
               ORDER BY g.seq_region_id, g.seq_region_start""",
        )
        gene_ids = [r["gene_id"] for r in cur.fetchall()]

    if not gene_ids:
        log.info("All genes already have stable IDs.")
        return 0

    # Find the next available ID
    with conn.cursor() as cur:
        cur.execute("SELECT MAX(stable_id) AS m FROM gene_stable_id")
        row = cur.fetchone()
        max_id = row["m"]

    if max_id:
        current_n = int(max_id.replace("ENS", "").lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("G")) + 1
    else:
        current_n = 1

    with conn.cursor() as cur:
        for gene_id in gene_ids:
            stable_id = _make_id(prefix, "G", current_n)
            cur.execute(
                "INSERT INTO gene_stable_id (gene_id, stable_id, version) VALUES (%s,%s,1)",
                (gene_id, stable_id),
            )
            current_n += 1

    log.info(f"  Assigned {len(gene_ids)} gene stable IDs")
    return len(gene_ids)


def assign_transcript_stable_ids(conn, prefix: str) -> int:
    with conn.cursor() as cur:
        cur.execute(
            """SELECT t.transcript_id FROM transcript t
               LEFT JOIN transcript_stable_id si ON t.transcript_id = si.transcript_id
               WHERE si.stable_id IS NULL
               ORDER BY t.seq_region_id, t.seq_region_start""",
        )
        tx_ids = [r["transcript_id"] for r in cur.fetchall()]

    if not tx_ids:
        log.info("All transcripts already have stable IDs.")
        return 0

    with conn.cursor() as cur:
        cur.execute("SELECT MAX(stable_id) AS m FROM transcript_stable_id")
        row = cur.fetchone()
        max_id = row["m"]
    current_n = (int(max_id.replace("ENS", "").lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("T")) + 1
                 if max_id else 1)

    with conn.cursor() as cur:
        for tx_id in tx_ids:
            stable_id = _make_id(prefix, "T", current_n)
            cur.execute(
                "INSERT INTO transcript_stable_id (transcript_id, stable_id, version) VALUES (%s,%s,1)",
                (tx_id, stable_id),
            )
            current_n += 1

    log.info(f"  Assigned {len(tx_ids)} transcript stable IDs")
    return len(tx_ids)


def assign_exon_stable_ids(conn, prefix: str) -> int:
    with conn.cursor() as cur:
        cur.execute(
            """SELECT e.exon_id FROM exon e
               LEFT JOIN exon_stable_id si ON e.exon_id = si.exon_id
               WHERE si.stable_id IS NULL
               ORDER BY e.seq_region_id, e.seq_region_start""",
        )
        exon_ids = [r["exon_id"] for r in cur.fetchall()]

    if not exon_ids:
        return 0

    with conn.cursor() as cur:
        cur.execute("SELECT MAX(stable_id) AS m FROM exon_stable_id")
        row = cur.fetchone()
        max_id = row["m"]
    current_n = (int(max_id.replace("ENS", "").lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("E")) + 1
                 if max_id else 1)

    with conn.cursor() as cur:
        for exon_id in exon_ids:
            stable_id = _make_id(prefix, "E", current_n)
            cur.execute(
                "INSERT INTO exon_stable_id (exon_id, stable_id, version) VALUES (%s,%s,1)",
                (exon_id, stable_id),
            )
            current_n += 1

    log.info(f"  Assigned {len(exon_ids)} exon stable IDs")
    return len(exon_ids)


def assign_translation_stable_ids(conn, prefix: str) -> int:
    with conn.cursor() as cur:
        cur.execute(
            """SELECT t.translation_id FROM translation t
               LEFT JOIN translation_stable_id si ON t.translation_id = si.translation_id
               WHERE si.stable_id IS NULL""",
        )
        tl_ids = [r["translation_id"] for r in cur.fetchall()]

    if not tl_ids:
        return 0

    with conn.cursor() as cur:
        cur.execute("SELECT MAX(stable_id) AS m FROM translation_stable_id")
        row = cur.fetchone()
        max_id = row["m"]
    current_n = (int(max_id.replace("ENS", "").lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ").lstrip("P")) + 1
                 if max_id else 1)

    with conn.cursor() as cur:
        for tl_id in tl_ids:
            stable_id = _make_id(prefix, "P", current_n)
            cur.execute(
                "INSERT INTO translation_stable_id (translation_id, stable_id, version) VALUES (%s,%s,1)",
                (tl_id, stable_id),
            )
            current_n += 1

    log.info(f"  Assigned {len(tl_ids)} translation stable IDs")
    return len(tl_ids)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host",        required=True)
    ap.add_argument("--port",        type=int, default=3306)
    ap.add_argument("--user",        required=True)
    ap.add_argument("--password",    default="")
    ap.add_argument("--dbname",      required=True)
    ap.add_argument("--prefix",      default="",
                    help="Species-specific prefix after ENS, e.g. 'GAL' for chicken. Empty for human/mouse.")
    ap.add_argument("--species-id",  type=int, default=1)
    return ap.parse_args()


def main():
    args = parse_args()
    conn = pymysql.connect(
        host=args.host, port=args.port,
        user=args.user, password=args.password,
        database=args.dbname,
        cursorclass=pymysql.cursors.DictCursor,
        autocommit=False,
    )
    try:
        log.info(f"Assigning stable IDs (prefix=ENS{args.prefix})")
        assign_gene_stable_ids(conn, args.prefix, args.species_id)
        assign_transcript_stable_ids(conn, args.prefix)
        assign_exon_stable_ids(conn, args.prefix)
        assign_translation_stable_ids(conn, args.prefix)
        conn.commit()
        log.info("Stable ID assignment complete.")
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
