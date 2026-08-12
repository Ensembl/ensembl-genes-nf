#!/usr/bin/env python3
"""Run PR 247's GFF feature loader against a sequence-free SQLite core stub.

This intentionally does not load FASTA or DNA. It creates only the coordinate
system, seq_region, and feature tables needed by load_gff_features_to_core,
then lets the upstream parser, insert code, and quality checks run unchanged.
It is a plumbing probe, not a production core builder.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
from dataclasses import replace
from pathlib import Path


def schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE TABLE coord_system (
            coord_system_id INTEGER PRIMARY KEY AUTOINCREMENT,
            species_id INTEGER NOT NULL DEFAULT 1,
            name TEXT NOT NULL,
            version TEXT,
            rank INTEGER,
            attrib TEXT
        );
        CREATE TABLE seq_region (
            seq_region_id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL,
            coord_system_id INTEGER NOT NULL,
            length INTEGER NOT NULL
        );
        CREATE TABLE analysis (
            analysis_id INTEGER PRIMARY KEY AUTOINCREMENT,
            logic_name TEXT,
            created TEXT,
            program TEXT
        );
        CREATE TABLE gene (
            gene_id INTEGER PRIMARY KEY,
            seq_region_id INTEGER,
            seq_region_start INTEGER,
            seq_region_end INTEGER,
            seq_region_strand INTEGER,
            biotype TEXT,
            analysis_id INTEGER,
            stable_id TEXT,
            display_xref_id INTEGER,
            source TEXT,
            canonical_transcript_id INTEGER
        );
        CREATE TABLE transcript (
            transcript_id INTEGER PRIMARY KEY,
            gene_id INTEGER,
            seq_region_id INTEGER,
            seq_region_start INTEGER,
            seq_region_end INTEGER,
            seq_region_strand INTEGER,
            biotype TEXT,
            analysis_id INTEGER,
            stable_id TEXT,
            display_xref_id INTEGER,
            source TEXT,
            canonical_translation_id INTEGER
        );
        CREATE TABLE exon (
            exon_id INTEGER PRIMARY KEY,
            seq_region_id INTEGER,
            seq_region_start INTEGER,
            seq_region_end INTEGER,
            seq_region_strand INTEGER,
            phase INTEGER,
            end_phase INTEGER,
            stable_id TEXT
        );
        CREATE TABLE exon_transcript (exon_id INTEGER, transcript_id INTEGER, rank INTEGER);
        CREATE TABLE translation (
            translation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            transcript_id INTEGER,
            start_exon_id INTEGER,
            end_exon_id INTEGER,
            seq_start INTEGER,
            seq_end INTEGER,
            stable_id TEXT
        );
        """
    )


def sequence_regions(connection: sqlite3.Connection, gffs: list[Path]) -> int:
    lengths: dict[str, int] = {}
    for gff in gffs:
        for raw in gff.open(encoding="utf-8"):
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) == 9:
                lengths[fields[0]] = max(lengths.get(fields[0], 0), int(fields[4]))
    connection.execute(
        "INSERT INTO coord_system(name, version, rank, attrib) VALUES (?, ?, ?, ?)",
        ("primary_assembly", "", 1, "default_version,sequence_level"),
    )
    coord_id = connection.execute(
        "SELECT coord_system_id FROM coord_system"
    ).fetchone()[0]
    connection.executemany(
        "INSERT INTO seq_region(name, coord_system_id, length) VALUES (?, ?, ?)",
        [(name, coord_id, length) for name, length in sorted(lengths.items())],
    )
    connection.commit()
    return len(lengths)


class CursorAdapter:
    def __init__(self, cursor: sqlite3.Cursor) -> None:
        self.cursor = cursor

    @property
    def lastrowid(self) -> int:
        return self.cursor.lastrowid

    def execute(self, statement: str, params=None):
        statement = statement.replace("%s", "?").replace("NOW()", "CURRENT_TIMESTAMP")
        return self.cursor.execute(statement, () if params is None else params)

    def executemany(self, statement: str, params):
        return self.cursor.executemany(statement.replace("%s", "?"), params)

    def fetchone(self):
        return self.cursor.fetchone()

    def fetchall(self):
        return self.cursor.fetchall()


class ConnectionAdapter:
    def __init__(self, connection: sqlite3.Connection) -> None:
        self.connection = connection

    def cursor(self) -> CursorAdapter:
        return CursorAdapter(self.connection.cursor())

    def commit(self) -> None:
        self.connection.commit()

    def rollback(self) -> None:
        self.connection.rollback()

    def close(self) -> None:
        # The upstream loader closes its DB connection before this probe prints
        # counts. Keep the SQLite handle open until the probe has reported them.
        pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", type=Path)
    parser.add_argument("--append-gff", type=Path, action="append", default=[])
    parser.add_argument("--pr247-root", type=Path, required=True)
    parser.add_argument("--source", default="ensembl")
    parser.add_argument("--append-source", action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(args.pr247_root / "src/python"))
    from ensembl.genes.ensembl_loading import gff_core_loader
    from ensembl.genes.ensembl_loading.gff_source_config import (
        ENSEMBL_GFF_CONFIG,
        get_source_config,
    )

    gff_paths = [args.gff, *args.append_gff]
    if args.output.exists():
        args.output.unlink()
    connection = sqlite3.connect(args.output)
    schema(connection)
    n_regions = sequence_regions(connection, gff_paths)
    gff_core_loader.connect_mysql = lambda **_: ConnectionAdapter(connection)
    def source_config(name: str):
        if name == "gencode50_projection":
            return replace(
            ENSEMBL_GFF_CONFIG,
            name="gencode50_projection",
            source_label="gencode50_projection",
            analysis_logic_name="gencode50_projection",
            analysis_program="GENCODE50_projection",
            gene_biotype_attribute="gene_type",
            transcript_biotype_attribute="transcript_type",
        )
        return get_source_config(name)

    source_configs = [source_config(args.source)]
    source_configs.extend(source_config(name) for name in args.append_source)
    if len(source_configs) != len(gff_paths):
        parser.error("--append-source must be supplied once for each --append-gff")

    summaries = []
    for gff_path, source_config in zip(gff_paths, source_configs):
        summaries.append(gff_core_loader.load_gff_features_to_core(
            gff_path=gff_path,
            db_name="sequence_free_probe",
            db_host="sqlite",
            db_user="probe",
            db_password="",
            db_port=0,
            coord_system_name="primary_assembly",
            source_config=source_config,
        ))
    print(f"seq_regions={n_regions}")
    print(summaries)
    for table in ("gene", "transcript", "exon", "exon_transcript", "translation"):
        count = connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"{table}={count}")
    connection.close()


if __name__ == "__main__":
    main()
