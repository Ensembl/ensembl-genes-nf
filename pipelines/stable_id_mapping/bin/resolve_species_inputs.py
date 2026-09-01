#!/usr/bin/env python3
"""Resolve all inputs for a single species stable-ID mapping run."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path

import pymysql

FTP_ROOT = Path("/nfs/ftp/public/ensemblorganisms")
GB_SERVER_HOST = os.environ["GBS1"]
GB_SERVER_PORT = int(os.environ["GBP1"])
ENSADMIN_PASSWORD = "ensembl"


class NoLiveReferenceError(Exception):
    pass


@dataclass
class SpeciesInputs:
    db_name: str
    target_chain: str
    target_version: int
    reference_chain: str
    reference_version: int
    binomial_name: str
    target_assembly_id: int
    reference_assembly_id: int
    target_genebuild_status_id: int
    ref_fasta: str
    ref_gff: str
    target_fasta: str
    target_gff: str
    mapping_session_id: int
    gene_range: str
    transcript_range: str
    translation_range: str


def _connect_ro():
    return pymysql.connect(
        host=GB_SERVER_HOST, port=GB_SERVER_PORT,
        user="ensro", password="",
        cursorclass=pymysql.cursors.DictCursor,
    )


def _connect_write():
    return pymysql.connect(
        host=GB_SERVER_HOST, port=GB_SERVER_PORT,
        user="ensadmin", password=ENSADMIN_PASSWORD,
        cursorclass=pymysql.cursors.DictCursor,
    )


def optional_path(value: str) -> Path | None:
    if value is None or value.strip() == "":
        return None
    return Path(value)


def parse_gca_accession(accession: str) -> tuple[str, int]:
    match = re.match(r"^(GCA_\d+)\.(\d+)$", accession)
    if not match:
        raise ValueError(f"Unrecognized GCA accession format: {accession}")
    return match.group(1), int(match.group(2))


def fetch_core_db_meta(db_name: str) -> dict[str, str]:
    conn = _connect_ro()
    try:
        with conn.cursor() as cur:
            cur.execute(f"USE `{db_name}`")
            cur.execute(
                "SELECT meta_key, meta_value FROM meta "
                "WHERE meta_key IN ('assembly.accession', 'species.scientific_name')"
            )
            rows = cur.fetchall()
    finally:
        conn.close()
    meta = {r["meta_key"]: r["meta_value"] for r in rows}
    for key in ("assembly.accession", "species.scientific_name"):
        if key not in meta:
            raise ValueError(f"Missing meta_key '{key}' in {db_name}")
    return meta


def resolve_reference_assembly(target_chain: str, target_version: int):
    conn = _connect_ro()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "USE gb_a_m_test"
            )
            cur.execute(
                "SELECT a.assembly_id, gs.genebuild_status_id "
                "FROM assembly a "
                "INNER JOIN genebuild_status gs ON gs.assembly_id = a.assembly_id "
                "WHERE a.gca_chain = %s AND a.gca_version = %s",
                (target_chain, target_version),
            )
            target_row = cur.fetchone()
            if target_row is None:
                raise ValueError(
                    f"No assembly/genebuild_status row for {target_chain}.{target_version}"
                )

            cur.execute(
                "SELECT a.assembly_id, a.gca_version "
                "FROM assembly a "
                "INNER JOIN genebuild_status gs ON gs.assembly_id = a.assembly_id "
                "WHERE a.gca_chain = %s AND a.gca_version < %s AND gs.gb_status = 'live' "
                "ORDER BY a.gca_version DESC LIMIT 1",
                (target_chain, target_version),
            )
            ref_row = cur.fetchone()
            if ref_row is None:
                raise NoLiveReferenceError(
                    f"No live reference version found for {target_chain} below version {target_version}"
                )
    finally:
        conn.close()

    return (
        target_chain,
        ref_row["gca_version"],
        ref_row["assembly_id"],
        target_row["assembly_id"],
        target_row["genebuild_status_id"],
    )


def resolve_stable_id_ranges(gca_accession: str) -> dict[str, str]:
    conn = _connect_ro()
    try:
        with conn.cursor() as cur:
            cur.execute("USE gb_a_m_test")
            cur.execute(
                "SELECT CONCAT(prefix.prefix, ':', stable.stable_space_start, "
                "'-', stable.stable_space_end) AS id_range "
                "FROM stable_space_species_log AS species_log "
                "INNER JOIN stable_space AS stable "
                "  ON stable.stable_space_id = species_log.stable_space_id "
                "INNER JOIN species_prefix AS prefix "
                "  ON prefix.lowest_taxon_id = species_log.lowest_taxon_id "
                "WHERE species_log.gca_accession = %s",
                (gca_accession,),
            )
            row = cur.fetchone()
    finally:
        conn.close()

    if row is None:
        raise ValueError(f"No stable-ID range found for {gca_accession}")

    prefix, span = row["id_range"].split(":")
    return {
        "gene": f"{prefix}G:{span}",
        "transcript": f"{prefix}T:{span}",
        "translation": f"{prefix}P:{span}",
    }


def insert_mapping_session(genebuild_status_id: int, reference_assembly_id: int) -> int:
    conn = _connect_write()
    try:
        with conn.cursor() as cur:
            cur.execute("USE gb_a_m_test")
            value = f"ref:{reference_assembly_id};date:{date.today().isoformat()}"
            cur.execute(
                "INSERT INTO annotation_events (genebuild_status_id, event, value) "
                "VALUES (%s, %s, %s)",
                (genebuild_status_id, "stable_id_mapping", value),
            )
            conn.commit()
            mapping_session_id = cur.lastrowid
    finally:
        conn.close()
    return mapping_session_id


def resolve_ftp_paths(binomial_name: str, gca_chain: str, gca_version: int) -> tuple[Path, Path]:
    formatted_name = binomial_name.replace(" ", "_")
    assembly_dir = FTP_ROOT / formatted_name / f"{gca_chain}.{gca_version}"
    genome_fasta = assembly_dir / "genome" / "unmasked.fa.gz"

    geneset_root = assembly_dir / "ensembl" / "geneset"
    if not geneset_root.is_dir():
        raise FileNotFoundError(f"No geneset directory found: {geneset_root}")
    latest_date_dir = sorted(geneset_root.iterdir())[-1]
    geneset_gff = latest_date_dir / "genes.gff3.gz"

    for path in (genome_fasta, geneset_gff):
        if not path.exists():
            raise FileNotFoundError(f"Expected file not found: {path}")

    return genome_fasta, geneset_gff


def resolve_species_inputs(
    db_name: str,
    target_fasta_override: Path | None = None,
    target_gff_override: Path | None = None,
    mapping_session_id_override: int | None = None,
    ref_fasta_override: Path | None = None,
    ref_gff_override: Path | None = None,
) -> SpeciesInputs:
    meta = fetch_core_db_meta(db_name)
    target_chain, target_version = parse_gca_accession(meta["assembly.accession"])
    binomial_name = meta["species.scientific_name"]

    (
        reference_chain, reference_version,
        reference_assembly_id, target_assembly_id, target_genebuild_status_id,
    ) = resolve_reference_assembly(target_chain, target_version)

    if target_fasta_override is not None and target_gff_override is not None:
        target_fasta, target_gff = target_fasta_override, target_gff_override
    elif target_fasta_override is None and target_gff_override is None:
        target_fasta, target_gff = resolve_ftp_paths(binomial_name, target_chain, target_version)
    else:
        raise ValueError("Provide both target_fasta and target_gff, or neither")
    
    if ref_fasta_override is not None and ref_gff_override is not None:
        ref_fasta, ref_gff = ref_fasta_override, ref_gff_override
    elif ref_fasta_override is None and ref_gff_override is None:
        ref_fasta, ref_gff = resolve_ftp_paths(binomial_name, reference_chain, reference_version)
    else:
        raise ValueError("Provide both ref_fasta and ref_gff, or neither")

    ranges = resolve_stable_id_ranges(f"{target_chain}.{target_version}")

    if mapping_session_id_override is not None:
        mapping_session_id = mapping_session_id_override
    else:
        mapping_session_id = insert_mapping_session(
            target_genebuild_status_id,
            reference_assembly_id,
        )

    return SpeciesInputs(
        db_name=db_name, target_chain=target_chain, target_version=target_version,
        reference_chain=reference_chain, reference_version=reference_version,
        binomial_name=binomial_name, target_assembly_id=target_assembly_id,
        reference_assembly_id=reference_assembly_id,
        target_genebuild_status_id=target_genebuild_status_id,
        ref_fasta=str(ref_fasta), ref_gff=str(ref_gff),
        target_fasta=str(target_fasta), target_gff=str(target_gff),
        mapping_session_id=mapping_session_id,
        gene_range=ranges["gene"], transcript_range=ranges["transcript"],
        translation_range=ranges["translation"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--target-fasta", default="")
    parser.add_argument("--target-gff", default="")
    parser.add_argument("--mapping-session-id", default="")
    parser.add_argument("--ref-fasta", default="")
    parser.add_argument("--ref-gff", default="")
    args = parser.parse_args()

    try:
        mapping_session_id = (
            int(args.mapping_session_id)
            if args.mapping_session_id.strip()
            else None
        )
        
        result = resolve_species_inputs(
            args.db_name,
            target_fasta_override=optional_path(args.target_fasta),
            target_gff_override=optional_path(args.target_gff),
            mapping_session_id_override=mapping_session_id,
            ref_fasta_override=optional_path(args.ref_fasta),
            ref_gff_override=optional_path(args.ref_gff),
        )
    except NoLiveReferenceError as exc:
        print(f"{exc} — no stable-ID mapping needed for {args.db_name}.")
        sys.exit(3)

    args.output_json.write_text(json.dumps(asdict(result), indent=2))


if __name__ == "__main__":
    main()
