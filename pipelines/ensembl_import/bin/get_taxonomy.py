"""
Fetch taxonomy info for a taxon_id from NCBI and build Ensembl-style
meta key/value pairs.
"""

from __future__ import annotations

import argparse
import logging
import os

from _mysq_helper import mysql_execute_query, mysql_fetch_data
from Bio import Entrez

logger = logging.getLogger(__name__)


def get_info_from_db(database: str, host: str, port: int, user: str, password: str) -> int:
    """
    Fetch a single taxonomy ID from the database.
    """
    queries = (
        """
        SELECT meta_value
        FROM meta
        WHERE meta_key = 'species.taxonomy_id'
        LIMIT 1
        """,
        """
        SELECT meta_value
        FROM meta
        WHERE meta_key = 'organism.taxonomy_id'
        LIMIT 1
        """,
    )
    for query in queries:
        results = mysql_fetch_data(query, database, host, port, user, password=password)
        if results:
            return int(results[0][0])

    raise ValueError("No taxonomy ID found in database")


def get_species_taxonomy_meta_from_ncbi(
    taxon_id: int,
    email: str,
    api_key: str | None,
) -> list[tuple[str, str]]:
    logger.info(f"Fetching taxonomy info from NCBI for taxon ID: {taxon_id}")

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    with Entrez.efetch(db="taxonomy", id=str(taxon_id), retmode="xml") as handle:
        records = Entrez.read(handle)

    if not records:
        raise ValueError(f"No taxonomy record found in NCBI for taxon_id {taxon_id}")

    rec = records[0]
    meta_pairs: list[tuple[str, str]] = [("species.taxonomy_id", str(taxon_id))]

    scientific_name = rec.get("ScientificName")
    if scientific_name:
        meta_pairs.append(("species.scientific_name", scientific_name))
        logger.info(f"species.scientific_name: {scientific_name}")

    common_name = rec.get("CommonName")
    if common_name:
        meta_pairs.append(("species.common_name", common_name))
        logger.info(f"species.common_name: {common_name}")

    lineage_ex = rec.get("LineageEx", [])
    reversed_lineage = list(reversed(lineage_ex))

    logger.info("Walking lineage")
    for ancestor in reversed_lineage:
        name = ancestor.get("ScientificName")
        rank = ancestor.get("Rank", "")

        if rank == "genus":
            continue

        # NCBI places this root node before Eukaryota in LineageEx, but it is
        # not part of the Ensembl species classification values.
        if name == "cellular organisms":
            continue

        if name:
            meta_pairs.append(("species.classification", name))
            logger.info(f"species.classification: {name}")

        if name == "Eukaryota" or rank == "superkingdom":
            logger.info("Eukaryota/superkingdom reached")
            break

    logger.info("Lineage walk complete")
    return meta_pairs


def load_to_db(
    database: str,
    host: str,
    port: int,
    user: str,
    password: str,
    meta_pairs: list[tuple[str, str]],
) -> None:
    """
    Load taxonomy data into the database.
    """
    logger.info("Loading to database")

    for key, value in meta_pairs:
        query = "INSERT IGNORE INTO meta (meta_key, meta_value, species_id) " "VALUES (%s, %s, 1)"
        mysql_execute_query(query, database, host, port, user, password=password, params=(key, value))
    logger.info("Load complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch taxonomy info from NCBI")
    parser.add_argument("--taxon_id", type=int, help="NCBI taxonomy ID")
    parser.add_argument("--email", type=str, help="Optional email address for NCBI Entrez")
    parser.add_argument(
        "--api_key",
        type=str,
        help="NCBI API key (optional). Defaults to environment variable NCBI_API_KEY",
    )
    parser.add_argument("--database", "--db_name", type=str, help="Database name")
    parser.add_argument("--host", "--db_host", type=str, help="Database host")
    parser.add_argument("--port", "--db_port", type=int, help="Database port")
    parser.add_argument("--user_w", "--db_user", type=str, help="Database user")
    parser.add_argument("--password", "--db_password", type=str, help="Database password")
    parser.add_argument("--read_user", "--db_read_user", type=str, help="Database read user")
    parser.add_argument(
        "--core_mode",
        action="store_true",
        help="Enable flag to run on ensembl core DB",
    )
    args = parser.parse_args()

    email = args.email or "dev@ensembl.org"
    api_key = args.api_key or os.environ.get("NCBI_API_KEY")

    if args.core_mode:
        logger.info("Running in core mode - getting taxon ID from database")
        if not args.database or not args.host or args.port is None or not args.read_user:
            raise ValueError("--database, --host, --port and --user are required in core mode")
        tax_id = get_info_from_db(args.database, args.host, args.port, args.read_user, password="")
        logger.info(f"Taxon ID found: {tax_id}")
    else:
        if args.taxon_id is None:
            raise ValueError("--taxon_id is required in standalone mode")
        tax_id = args.taxon_id
        logger.info(f"Running in standalone mode for taxon ID: {tax_id}")

    meta = get_species_taxonomy_meta_from_ncbi(tax_id, email, api_key)

    for key, value in meta:
        print(key, "=>", value)

    if args.core_mode:
        logger.info("Running in core mode - loading to database")
        load_to_db(args.database, args.host, args.port, args.user_w, args.password, meta)
