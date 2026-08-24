from __future__ import annotations

import argparse
import logging
import re
from typing import TextIO

import pymysql

logger = logging.getLogger(__name__)

# We store thes PEPSTATS codes in the core

PEPSTATS_CODES = [
    "NumResidues",
    "MolecularWeight",
    "AvgResWeight",
    "Charge",
    "IsoPoint",
]


def parse_pepstats_stream(handle: TextIO) -> dict[str, dict[str, str]]:
    """
    Parse EMBOSS pepstats output from a text stream.

    Returns:
        {
          "<translation_id>": {
              "NumResidues": "<int as string>",
              "MolecularWeight": "<float as string>",
              "AvgResWeight": "<float as string>",
              "Charge": "<float as string>",
              "IsoPoint": "<float as string>",
          },
          ...
        }
    """
    attribs: dict[str, dict[str, str]] = {}
    tid: str | None = None

    # Precompile regexes equivalent to the Perl ones
    re_tid = re.compile(r"PEPSTATS of ([^ ]+)")
    re_mw_res = re.compile(r"^Molecular weight = (\S+)\s+Residues = (\d+).*")
    re_avg_charge = re.compile(r"^Average(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)")
    re_iso = re.compile(r"^Isoelectric(\s+)(\S+)(\s+)=(\s+)(\S+)")

    for line in handle:
        line = line.rstrip("\n")

        m_tid = re_tid.search(line)
        if m_tid:
            tid = m_tid.group(1)
            if tid not in attribs:
                attribs[tid] = {}
            continue

        if tid is None:
            # Ignore lines until we see "PEPSTATS of ..."
            continue

        m_mw = re_mw_res.match(line)
        if m_mw:
            mw = m_mw.group(1)
            residues = m_mw.group(2)
            attribs.setdefault(tid, {})
            attribs[tid]["NumResidues"] = residues
            attribs[tid]["MolecularWeight"] = mw
            continue

        m_avg = re_avg_charge.match(line)
        if m_avg:
            # Perl groups:
            # 1: spaces, 2: word, 3: spaces, 4: word, 5: spaces,
            # 6: spaces, 7: value1, 8: spaces, 9: word, 10: spaces,
            # 11: spaces, 12: value2
            avg_res_weight = m_avg.group(7)
            charge = m_avg.group(12)
            attribs.setdefault(tid, {})
            attribs[tid]["AvgResWeight"] = avg_res_weight
            attribs[tid]["Charge"] = charge
            continue

        m_iso = re_iso.match(line)
        if m_iso:
            iso_point = m_iso.group(5)
            attribs.setdefault(tid, {})
            attribs[tid]["IsoPoint"] = iso_point
            continue
        logger.debug(f"Line not matched: {line}")

    logger.info(f"Parsed {len(attribs)} attribs from pepstats file: {handle}")
    return attribs


def parse_pepstats_file(path: str) -> dict[str, dict[str, str]]:
    """
    Convenience wrapper: parse pepstats output from a file path.
    """
    with open(path, "r", encoding="utf-8") as fh:
        logger.info(f"Parsing pepstats file: {path}")
        return parse_pepstats_stream(fh)


def get_attrib_type_ids(conn) -> dict[str, int]:
    """
    Fetch attrib_type_id for the five pepstats codes from the core DB.

    Returns:
        { "NumResidues": 123, "MolecularWeight": 124, ... }
    """

    code_to_id: dict[str, int] = {}

    sql = f"""
        SELECT code, attrib_type_id
        FROM attrib_type
        WHERE code IN ({",".join(["%s"] * len(PEPSTATS_CODES))})
    """

    with conn.cursor() as cur:
        cur.execute(sql, PEPSTATS_CODES)
        for code, attrib_type_id in cur.fetchall():
            code_to_id[code] = attrib_type_id
    logger.info(f"Fetched {len(code_to_id)} attrib_type_id rows from core DB")

    missing = [c for c in PEPSTATS_CODES if c not in code_to_id]
    if missing:
        raise RuntimeError(f"Missing attrib_type rows for pepstats codes: {', '.join(missing)}")

    return code_to_id


def get_translation_ids(conn, stable_ids: list[str]) -> dict[str, int]:
    """Resolve FASTA translation stable IDs to core translation IDs."""
    if not stable_ids:
        return {}

    placeholders = ",".join(["%s"] * len(stable_ids))
    sql = f"""
        SELECT stable_id, translation_id
        FROM translation
        WHERE stable_id IN ({placeholders})
    """
    with conn.cursor() as cur:
        cur.execute(sql, stable_ids)
        return {stable_id: translation_id for stable_id, translation_id in cur.fetchall()}


def apply_pepstats_to_core(
    conn,
    pepstats_results: dict[str, dict[str, str]],
    species_id: int,
) -> None:
    """
    Apply parsed pepstats results to a core DB using pymysql.

    Args:
        conn: open pymysql connection to the *core* DB.
        pepstats_results:
            {
              "<translation_id>": {
                  "NumResidues": "317",
                  "MolecularWeight": "35513.44",
                  "AvgResWeight": "112.0",
                  "Charge": "-3.2",
                  "IsoPoint": "6.45",
              },
              ...
            }
        species_id: coord_system.species_id to restrict deletions.
    """

    code_to_id = get_attrib_type_ids(conn)
    stable_ids = [tid for tid in pepstats_results if not tid.isdigit()]
    translation_ids = get_translation_ids(conn, stable_ids)

    # 1. Delete existing pepstats attributes for this species
    delete_sql = f"""
        DELETE ta
        FROM translation_attrib ta
        JOIN attrib_type at  ON at.attrib_type_id = ta.attrib_type_id
        JOIN translation tl  ON tl.translation_id = ta.translation_id
        JOIN transcript tr   ON tr.transcript_id = tl.transcript_id
        JOIN seq_region s    ON s.seq_region_id = tr.seq_region_id
        JOIN coord_system c  ON c.coord_system_id = s.coord_system_id
        WHERE c.species_id = %s
          AND at.code IN ({",".join(["%s"] * len(PEPSTATS_CODES))})
    """

    with conn.cursor() as cur:
        cur.execute(delete_sql, [species_id, *PEPSTATS_CODES])
        logger.info(f"Cleaned {cur.rowcount} pepstats rows for species_id: {species_id}")

    # 2. Insert new rows into translation_attrib
    insert_sql = """
        INSERT INTO translation_attrib (translation_id, attrib_type_id, value)
        VALUES (%s, %s, %s)
    """

    data = []
    for tid_str, metrics in pepstats_results.items():
        if tid_str.isdigit():
            translation_id = int(tid_str)
        else:
            translation_id = translation_ids.get(tid_str)
            if translation_id is None:
                logger.warning("Translation stable ID not found in core: %s", tid_str)
                continue
        for code, value in metrics.items():
            if code not in code_to_id:
                continue
            attrib_type_id = code_to_id[code]
            data.append((translation_id, attrib_type_id, str(value)))
    logger.info(f"Prepared {len(data)} pepstats rows for insertion")

    if data:
        with conn.cursor() as cur:
            cur.executemany(insert_sql, data)
            logger.info(f"Inserted {cur.rowcount} pepstats rows into translation_attrib")

    conn.commit()
    logger.info("Database commit completed")


def main():
    """
    Parse a Pepstats file and add the data to the core database.
    """
    parser = argparse.ArgumentParser(description="Add Pepstats to core database")
    parser.add_argument("--pepstats_file", type=str, help="Path to the Pepstats file")
    parser.add_argument(
        "--species_id", default=1, type=int, help="Species ID to associate with the Pepstats, default is 1"
    )
    parser.add_argument("--db_host", type=str, help="Core database host")
    parser.add_argument("--db_port", type=int, help="Core database port")
    parser.add_argument("--db_user", type=str, help="Core database user (write access)")
    parser.add_argument("--db_password", type=str, help="Core database password (write access)")
    parser.add_argument("--db_name", type=str, help="Core database name")

    args = parser.parse_args()
    logger.info("Start parsing pepstats script")
    logger.info(f"Using predefined set {PEPSTATS_CODES}")

    parsed_pepstats = parse_pepstats_file(args.pepstats_file)
    conn = pymysql.connect(
        host=args.db_host,
        port=args.db_port,
        user=args.db_user,
        password=args.db_password,
        database=args.db_name,
    )
    logger.info(f"Using core database: {args.db_name} on {args.db_host}")

    apply_pepstats_to_core(conn, parsed_pepstats, args.species_id)


if __name__ == "__main__":
    main()
