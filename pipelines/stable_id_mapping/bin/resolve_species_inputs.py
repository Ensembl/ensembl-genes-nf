#!/usr/bin/env python3
"""Resolve all inputs for a single species stable-ID mapping run."""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass
from datetime import date, datetime
from pathlib import Path

import pymysql

from stable_id_mapping.ids import (
    exon_range_from_gene_range,
    parse_id_range,
)
from stable_id_mapping.range_check import (
    DuplicateStableIdError,
    check_stable_id_population,
    find_duplicate_stable_ids,
    load_current_features,
)


FTP_ROOT = Path("/nfs/ftp/public/ensemblorganisms")
PRE_RELEASE_ROOT = Path("/nfs/ftp/public/databases/ensembl/pre-release")
REG_SERVER_HOST: str
REG_SERVER_PORT: int
GB_SERVER_HOST: str
GB_SERVER_PORT: int
ENSADMIN_PASSWORD: str


class NoLiveReferenceError(Exception):
    pass


@dataclass
class SpeciesInputs:
    db_name: str
    requested_mode: str
    effective_mode: str
    target_chain: str
    target_version: int
    reference_chain: str | None
    reference_version: int | None
    binomial_name: str
    target_assembly_id: int
    reference_assembly_id: int | None
    target_genebuild_status_id: int
    ref_fasta: str | None
    ref_gff: str | None
    target_fasta: str | None
    target_gff: str | None
    mapping_session_id: int | None
    gene_range: str
    transcript_range: str
    translation_range: str


def _connect_reg_ro():
    return pymysql.connect(
        host=REG_SERVER_HOST, port=REG_SERVER_PORT,
        user="ensro", password="",
        cursorclass=pymysql.cursors.DictCursor,
    )


def _connect_reg_write():
    return pymysql.connect(
        host=REG_SERVER_HOST, port=REG_SERVER_PORT,
        user="ensadmin", password=ENSADMIN_PASSWORD,
        cursorclass=pymysql.cursors.DictCursor,
    )


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


def resolve_target_assembly(
    target_chain: str,
    target_version: int,
    assembly_metadata_db: str,
) -> tuple[int, int]:
    conn = _connect_reg_ro()
    try:
        with conn.cursor() as cur:
            cur.execute(f"USE {assembly_metadata_db}")
            cur.execute(
                "SELECT a.assembly_id, gs.genebuild_status_id "
                "FROM assembly a "
                "INNER JOIN genebuild_status gs "
                "  ON gs.assembly_id = a.assembly_id "
                "WHERE a.gca_chain = %s AND a.gca_version = %s",
                (target_chain, target_version),
            )
            row = cur.fetchone()
    finally:
        conn.close()

    if row is None:
        raise ValueError(
            f"No assembly/genebuild_status row for "
            f"{target_chain}.{target_version}"
        )

    return row["assembly_id"], row["genebuild_status_id"]


def resolve_live_reference(
    target_chain: str,
    target_version: int,
    assembly_metadata_db: str,
) -> tuple[str, int, int] | None:
    conn = _connect_reg_ro()
    try:
        with conn.cursor() as cur:
            cur.execute(f"USE {assembly_metadata_db}")
            cur.execute(
                "SELECT a.assembly_id, a.gca_version "
                "FROM assembly a "
                "INNER JOIN genebuild_status gs "
                "  ON gs.assembly_id = a.assembly_id "
                "WHERE a.gca_chain = %s "
                "  AND a.gca_version < %s "
                "  AND gs.gb_status = 'live' "
                "ORDER BY a.gca_version DESC "
                "LIMIT 1",
                (target_chain, target_version),
            )
            row = cur.fetchone()
    finally:
        conn.close()

    if row is None:
        return None

    return (
        target_chain,
        row["gca_version"],
        row["assembly_id"],
    )


def resolve_stable_id_ranges(
    gca_accession: str,
    assembly_metadata_db: str,
) -> dict[str, str]:
    conn = _connect_reg_ro()
    try:
        with conn.cursor() as cur:
            cur.execute(f"USE {assembly_metadata_db}")
            cur.execute(
                "SELECT prefix.prefix AS prefix, "
                "stable.stable_space_start AS range_start, "
                "stable.stable_space_end AS range_end "
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

    start = int(row["range_start"])
    end = int(row["range_end"])
    maximum = (10**11) - 1

    if start < 0 or end < start or end > maximum:
        raise ValueError(
            f"Invalid 11-digit stable-ID range for {gca_accession}: "
            f"{start}-{end}"
        )

    prefix = row["prefix"]
    span = f"{start:011d}-{end:011d}"

    return {
        "gene": f"{prefix}G:{span}",
        "transcript": f"{prefix}T:{span}",
        "translation": f"{prefix}P:{span}",
    }


def stable_ids_require_reassignment(
    db_name: str,
    ranges: dict[str, str],
) -> bool:
    gene_range = parse_id_range(ranges["gene"])

    expected_ranges = {
        "gene": gene_range,
        "transcript": parse_id_range(ranges["transcript"]),
        "translation": parse_id_range(ranges["translation"]),
        "exon": exon_range_from_gene_range(gene_range),
    }

    populations = load_current_features(db_name)

    all_stable_ids = [
        row["stable_id"]
        for rows in populations.values()
        for row in rows
    ]
    duplicates = find_duplicate_stable_ids(all_stable_ids)

    if duplicates:
        examples = ", ".join(duplicates[:20])
        raise DuplicateStableIdError(
            "Duplicate stable IDs detected across core feature tables "
            f"({len(duplicates)} duplicated IDs). Examples: {examples}"
        )

    checks = {
        feature_type: check_stable_id_population(
            (
                row["stable_id"]
                for row in populations[feature_type]
            ),
            expected_range,
        )
        for feature_type, expected_range in expected_ranges.items()
    }

    agreeing = sum(check.agreeing for check in checks.values())
    disagreeing = sum(check.disagreeing for check in checks.values())

    if agreeing and disagreeing:
        print(
            "\n*** WARNING: MIXED STABLE-ID POPULATION DETECTED ***\n"
            f"Database: {db_name}\n"
            f"IDs agreeing with registry: {agreeing}\n"
            f"IDs requiring reassignment: {disagreeing}\n"
            "The complete reassignment SQL will be generated for review.",
        )

    return disagreeing > 0


def insert_mapping_session(
    genebuild_status_id: int,
    reference_assembly_id: int,
    assembly_metadata_db: str,
    gca_chain: str,
) -> int:
    conn = _connect_reg_write()
    try:
        with conn.cursor() as cur:
            cur.execute(f"USE {assembly_metadata_db}")
            value = f"ref:{reference_assembly_id}({gca_chain});date:{datetime.now().isoformat(timespec='minutes')};not_executed"
            cur.execute(
                "INSERT ignore INTO annotation_events (genebuild_status_id, event, value) "
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


def resolve_pre_release_paths(
    binomial_name: str,
    gca_chain: str,
    gca_version: int,
) -> tuple[Path, Path]:
    formatted_name = binomial_name.replace(" ", "_")
    assembly_dir = (
        PRE_RELEASE_ROOT
        / formatted_name
        / f"{gca_chain}.{gca_version}"
    )

    if not assembly_dir.is_dir():
        raise FileNotFoundError(
            f"No pre-release assembly directory found: {assembly_dir}"
        )

    fasta_matches = sorted(
        assembly_dir.glob("*.dna.softmasked.fa.gz")
    )
    gff_matches = sorted(
        assembly_dir.glob("*.gff3.gz")
    )

    if len(fasta_matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one pre-release softmasked FASTA in "
            f"{assembly_dir}, found {len(fasta_matches)}"
        )

    if len(gff_matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one pre-release GFF3 in "
            f"{assembly_dir}, found {len(gff_matches)}"
        )

    return fasta_matches[0], gff_matches[0]


def resolve_species_inputs(
    db_name: str,
    requested_mode: str = "auto",
    assembly_metadata_db: str = "gb_assembly_metadata",
    target_fasta_override: Path | None = None,
    target_gff_override: Path | None = None,
    mapping_session_id_override: int | None = None,
    ref_fasta_override: Path | None = None,
    ref_gff_override: Path | None = None,
) -> SpeciesInputs:
    if requested_mode not in {"auto", "map", "reassign"}:
        raise ValueError(
            f"Unknown mode {requested_mode!r}; "
            "expected auto, map, or reassign"
        )

    meta = fetch_core_db_meta(db_name)
    target_chain, target_version = parse_gca_accession(
        meta["assembly.accession"]
    )
    binomial_name = meta["species.scientific_name"]

    (
        target_assembly_id,
        target_genebuild_status_id,
    ) = resolve_target_assembly(
        target_chain,
        target_version,
        assembly_metadata_db,
    )

    ranges = resolve_stable_id_ranges(
        f"{target_chain}.{target_version}",
        assembly_metadata_db
    )


    reference = None

    if requested_mode in {"auto", "map"}:
        reference = resolve_live_reference(
            target_chain,
            target_version,
            assembly_metadata_db,
        )

    if requested_mode == "map":
        if reference is None:
            raise NoLiveReferenceError(
                f"Mode 'map' requested, but no live reference version "
                f"exists for {target_chain} below version {target_version}"
            )
        effective_mode = "map"
    
    elif requested_mode == "reassign":
        stable_ids_require_reassignment(db_name, ranges)
        effective_mode = "reassign"
    
    elif reference is not None:
        effective_mode = "map"
    
    elif stable_ids_require_reassignment(db_name, ranges):
        effective_mode = "reassign"
    
    else:
        effective_mode = "no_action"

    if effective_mode == "map":
        assert reference is not None
    
        (
            reference_chain,
            reference_version,
            reference_assembly_id,
        ) = reference
    
        if target_fasta_override is not None and target_gff_override is not None:
            target_fasta = target_fasta_override
            target_gff = target_gff_override
        elif target_fasta_override is None and target_gff_override is None:
            try:
                target_fasta, target_gff = resolve_ftp_paths(
                    binomial_name,
                    target_chain,
                    target_version,
                )
            except FileNotFoundError as standard_error:
                try:
                    (
                        target_fasta,
                        target_gff,
                    ) = resolve_pre_release_paths(
                        binomial_name,
                        target_chain,
                        target_version,
                    )
                except FileNotFoundError as pre_release_error:
                    raise FileNotFoundError(
                        "Target files were unavailable in both the "
                        "standard and pre-release locations.\n"
                        f"Standard lookup: {standard_error}\n"
                        f"Pre-release lookup: {pre_release_error}"
                    ) from pre_release_error

                print(
                    "Target files were not found in the standard FTP "
                    f"tree; using pre-release files from "
                    f"{target_fasta.parent}"
                )
        else:
            raise ValueError(
                "Provide both target_fasta and target_gff, or neither"
            )
    
        if ref_fasta_override is not None and ref_gff_override is not None:
            ref_fasta = ref_fasta_override
            ref_gff = ref_gff_override
        elif ref_fasta_override is None and ref_gff_override is None:
            ref_fasta, ref_gff = resolve_ftp_paths(
                binomial_name,
                reference_chain,
                reference_version,
            )
        else:
            raise ValueError(
                "Provide both ref_fasta and ref_gff, or neither"
            )
    
        if mapping_session_id_override is not None:
            mapping_session_id = mapping_session_id_override
        else:
            mapping_session_id = insert_mapping_session(
                target_genebuild_status_id,
                reference_assembly_id,
                assembly_metadata_db,
                target_chain,
            )
    else:
        reference_chain = None
        reference_version = None
        reference_assembly_id = None
        target_fasta = None
        target_gff = None
        ref_fasta = None
        ref_gff = None
        mapping_session_id = None

    return SpeciesInputs(
        db_name=db_name,
        requested_mode=requested_mode,
        effective_mode=effective_mode,
        target_chain=target_chain,
        target_version=target_version,
        reference_chain=reference_chain,
        reference_version=reference_version,
        binomial_name=binomial_name,
        target_assembly_id=target_assembly_id,
        reference_assembly_id=reference_assembly_id,
        target_genebuild_status_id=target_genebuild_status_id,
        ref_fasta=str(ref_fasta) if ref_fasta is not None else None,
        ref_gff=str(ref_gff) if ref_gff is not None else None,
        target_fasta=(
            str(target_fasta)
            if target_fasta is not None
            else None
        ),
        target_gff=(
            str(target_gff)
            if target_gff is not None
            else None
        ),
        mapping_session_id=mapping_session_id,
        gene_range=ranges["gene"],
        transcript_range=ranges["transcript"],
        translation_range=ranges["translation"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--assembly-metadata-db", required=True)
    parser.add_argument(
        "--mode",
        choices=("auto", "map", "reassign"),
        default="auto",
    )
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--target-fasta", default="")
    parser.add_argument("--target-gff", default="")
    parser.add_argument("--mapping-session-id", default="")
    parser.add_argument("--ref-fasta", default="")
    parser.add_argument("--ref-gff", default="")
    parser.add_argument("--reg-host", required=True)
    parser.add_argument("--reg-port", required=True, type=int)
    parser.add_argument("--gb-host", required=True)
    parser.add_argument("--gb-port", required=True, type=int)
    parser.add_argument("--ensadmin-password", required=True)
    args = parser.parse_args()

    global REG_SERVER_HOST, REG_SERVER_PORT, GB_SERVER_HOST, GB_SERVER_PORT, ENSADMIN_PASSWORD
    REG_SERVER_HOST = args.reg_host
    REG_SERVER_PORT = args.reg_port
    GB_SERVER_HOST = args.gb_host
    GB_SERVER_PORT = args.gb_port
    ENSADMIN_PASSWORD = args.ensadmin_password
    
    try:
        mapping_session_id = (
            int(args.mapping_session_id)
            if args.mapping_session_id.strip()
            else None
        )
        
        result = resolve_species_inputs(
            args.db_name,
            requested_mode=args.mode,
            target_fasta_override=optional_path(args.target_fasta),
            target_gff_override=optional_path(args.target_gff),
            mapping_session_id_override=mapping_session_id,
            ref_fasta_override=optional_path(args.ref_fasta),
            ref_gff_override=optional_path(args.ref_gff),
        )
    except NoLiveReferenceError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    args.output_json.write_text(json.dumps(asdict(result), indent=2))


if __name__ == "__main__":
    main()
