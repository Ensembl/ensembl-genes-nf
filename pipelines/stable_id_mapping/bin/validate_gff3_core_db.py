#!/usr/bin/env python3
"""Compare stable IDs, versions, locations, and parent links in a GFF3 and core DB."""

from __future__ import annotations

import argparse
import getpass
import gzip
import sys
from pathlib import Path

import pymysql


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a rewritten target GFF3 against a core database."
    )
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--host", required=True)
    parser.add_argument("--port", type=int, default=3306)
    parser.add_argument("--user", default="ensadmin")
    parser.add_argument(
        "--password",
        help="MySQL password. If omitted, it is requested interactively.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_attributes(text: str) -> dict[str, str]:
    attributes = {}

    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Malformed GFF3 attribute: {item!r}")
        key, value = item.split("=", 1)
        attributes[key] = value

    return attributes


def bare_id(value: str) -> str:
    """Remove an optional GFF3 namespace, such as gene: or transcript:."""
    return value.split(":", 1)[-1]


def add_record(
    records: dict[str, dict[str, str]],
    stable_id: str,
    record: dict[str, str],
    line_number: int,
) -> None:
    previous = records.get(stable_id)

    if previous is None:
        records[stable_id] = record
    elif previous != record:
        raise ValueError(
            f"line {line_number}: inconsistent repeated record for {stable_id}"
        )


def load_gff_records(
    path: Path,
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
]:
    genes: dict[str, dict[str, str]] = {}
    transcripts: dict[str, dict[str, str]] = {}
    translations: dict[str, dict[str, str]] = {}

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"line {line_number}: expected nine GFF3 columns"
                )

            seqid, _, _, start, end, _, strand, _, attribute_text = fields
            attributes = parse_attributes(attribute_text)

            location = {
                "version": attributes.get("version", ""),
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
            }

            if "gene_id" in attributes:
                stable_id = bare_id(attributes["gene_id"])
                add_record(genes, stable_id, location, line_number)

            elif "transcript_id" in attributes:
                stable_id = bare_id(attributes["transcript_id"])
                parent = attributes.get("Parent", "")
                record = location | {
                    "gene_id": bare_id(parent),
                }
                add_record(transcripts, stable_id, record, line_number)

            elif "protein_id" in attributes:
                stable_id = bare_id(attributes["protein_id"])
                parent = attributes.get("Parent", "")
                record = {
                    "version": attributes.get("version", ""),
                    "transcript_id": bare_id(parent),
                }
                add_record(translations, stable_id, record, line_number)

    return genes, transcripts, translations


def database_strand_to_gff(value: int) -> str:
    strands = {
        1: "+",
        -1: "-",
        0: ".",
    }
    return strands[int(value)]


def load_database_records(
    connection: pymysql.Connection,
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
]:
    queries = {
        "genes": """
            SELECT
                gene.stable_id,
                gene.version,
                seq_region.name,
                gene.seq_region_start,
                gene.seq_region_end,
                gene.seq_region_strand
            FROM gene
            JOIN seq_region USING (seq_region_id)
            WHERE gene.stable_id IS NOT NULL
        """,
        "transcripts": """
            SELECT
                transcript.stable_id,
                transcript.version,
                seq_region.name,
                transcript.seq_region_start,
                transcript.seq_region_end,
                transcript.seq_region_strand,
                gene.stable_id
            FROM transcript
            JOIN gene USING (gene_id)
            JOIN seq_region
                ON transcript.seq_region_id = seq_region.seq_region_id
            WHERE transcript.stable_id IS NOT NULL
        """,
        "translations": """
            SELECT
                translation.stable_id,
                translation.version,
                transcript.stable_id
            FROM translation
            JOIN transcript USING (transcript_id)
            WHERE translation.stable_id IS NOT NULL
        """,
    }

    database_records = {}

    with connection.cursor() as cursor:
        for record_type, query in queries.items():
            cursor.execute(query)
            rows = cursor.fetchall()
            records = {}

            for row in rows:
                if record_type == "genes":
                    stable_id, version, seqid, start, end, strand = row
                    records[stable_id] = {
                        "version": str(version),
                        "seqid": seqid,
                        "start": str(start),
                        "end": str(end),
                        "strand": database_strand_to_gff(strand),
                    }

                elif record_type == "transcripts":
                    (
                        stable_id,
                        version,
                        seqid,
                        start,
                        end,
                        strand,
                        gene_id,
                    ) = row
                    records[stable_id] = {
                        "version": str(version),
                        "seqid": seqid,
                        "start": str(start),
                        "end": str(end),
                        "strand": database_strand_to_gff(strand),
                        "gene_id": gene_id,
                    }

                else:
                    stable_id, version, transcript_id = row
                    records[stable_id] = {
                        "version": str(version),
                        "transcript_id": transcript_id,
                    }

            database_records[record_type] = records

    return (
        database_records["genes"],
        database_records["transcripts"],
        database_records["translations"],
    )


def compare_records(
    label: str,
    gff_records: dict[str, dict[str, str]],
    database_records: dict[str, dict[str, str]],
) -> list[str]:
    errors = []
    gff_ids = set(gff_records)
    database_ids = set(database_records)

    only_in_gff = gff_ids - database_ids
    only_in_database = database_ids - gff_ids

    if only_in_gff:
        errors.append(
            f"{label}: {len(only_in_gff)} ID(s) exist only in GFF3 "
            f"(examples: {', '.join(sorted(only_in_gff)[:5])})"
        )

    if only_in_database:
        errors.append(
            f"{label}: {len(only_in_database)} ID(s) exist only in database "
            f"(examples: {', '.join(sorted(only_in_database)[:5])})"
        )

    for stable_id in sorted(gff_ids & database_ids):
        if gff_records[stable_id] != database_records[stable_id]:
            errors.append(
                f"{label}: {stable_id} differs\n"
                f"    GFF3:     {gff_records[stable_id]}\n"
                f"    database: {database_records[stable_id]}"
            )

    return errors


def main() -> int:
    args = parse_args()
    password = args.password or getpass.getpass("MySQL password: ")

    gff_genes, gff_transcripts, gff_translations = load_gff_records(args.gff)

    connection = pymysql.connect(
        host=args.host,
        port=args.port,
        user=args.user,
        password=password,
        database=args.db_name,
    )

    try:
        db_genes, db_transcripts, db_translations = load_database_records(
            connection
        )
    finally:
        connection.close()

    print("GFF3 / CORE DATABASE VALIDATION")
    print(f"GFF3:     {args.gff}")
    print(f"Database: {args.db_name}")
    print(f"Genes:        GFF3={len(gff_genes)} DB={len(db_genes)}")
    print(
        f"Transcripts:  GFF3={len(gff_transcripts)} "
        f"DB={len(db_transcripts)}"
    )
    print(
        f"Translations: GFF3={len(gff_translations)} "
        f"DB={len(db_translations)}"
    )

    errors = []
    errors.extend(compare_records("genes", gff_genes, db_genes))
    errors.extend(compare_records("transcripts", gff_transcripts, db_transcripts))
    errors.extend(compare_records("translations", gff_translations, db_translations))

    if errors:
        print(f"\nFAILED: {len(errors)} difference(s)", file=sys.stderr)
        for error in errors[:20]:
            print(f"\n{error}", file=sys.stderr)
        if len(errors) > 20:
            print(
                f"\n... and {len(errors) - 20} more",
                file=sys.stderr,
            )
        return 1

    print("\nPASSED: GFF3 and core database agree for all stable IDs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())