#!/usr/bin/env python3
"""Select and store the sample location metadata for a core database.

This follows the selection rules used by Ensembl's
HiveAddPlaceholderLocation runnable.  The database is expected to have been
created and loaded before this script is run.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from typing import Optional, Sequence

import pymysql

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class SampleGene:
    gene_stable_id: str
    transcript_stable_id: str
    seq_region_name: str
    gene_start: int
    gene_end: int

    @property
    def location(self) -> str:
        return f"{self.seq_region_name}:{self.gene_start}-{self.gene_end}"


def longest_seq_regions(connection: pymysql.connections.Connection) -> Sequence[int]:
    """Return the IDs of the ten longest sequence regions."""
    with connection.cursor() as cursor:
        cursor.execute("SELECT seq_region_id FROM seq_region ORDER BY length DESC LIMIT 10")
        return [row[0] for row in cursor.fetchall()]


def select_supported_transcript(
    connection: pymysql.connections.Connection, seq_region_id: int
) -> Optional[SampleGene]:
    """Select the first supported protein-coding transcript on a region."""
    query = """
        SELECT g.stable_id, t.stable_id, sr.name,
               g.seq_region_start, g.seq_region_end
          FROM transcript AS t
          JOIN gene AS g ON g.gene_id = t.gene_id
          JOIN seq_region AS sr ON sr.seq_region_id = g.seq_region_id
         WHERE t.seq_region_id = %s
           AND t.biotype = 'protein_coding'
           AND EXISTS (
               SELECT 1
                 FROM transcript_supporting_feature AS tsf
                WHERE tsf.transcript_id = t.transcript_id
                  AND (
                      (tsf.feature_type = 'dna_align_feature' AND EXISTS (
                          SELECT 1 FROM dna_align_feature AS daf
                           WHERE daf.dna_align_feature_id = tsf.feature_id
                             AND daf.hcoverage >= 99
                             AND daf.perc_ident >= 75
                      ))
                      OR
                      (tsf.feature_type = 'protein_align_feature' AND EXISTS (
                          SELECT 1 FROM protein_align_feature AS paf
                           WHERE paf.protein_align_feature_id = tsf.feature_id
                             AND paf.hcoverage >= 99
                             AND paf.perc_ident >= 75
                      ))
                  )
           )
         ORDER BY t.transcript_id
         LIMIT 1
    """
    with connection.cursor() as cursor:
        cursor.execute(query, (seq_region_id,))
        row = cursor.fetchone()
    return SampleGene(*row) if row else None


def select_largest_gene(
    connection: pymysql.connections.Connection, seq_region_id: int
) -> Optional[SampleGene]:
    """Select the protein-coding gene with the most exons on a region."""
    query = """
        SELECT g.stable_id, t.stable_id, sr.name,
               g.seq_region_start, g.seq_region_end
          FROM gene AS g
          JOIN transcript AS t ON t.transcript_id = g.canonical_transcript_id
          JOIN seq_region AS sr ON sr.seq_region_id = g.seq_region_id
          LEFT JOIN exon_transcript AS et ON et.transcript_id = t.transcript_id
         WHERE g.seq_region_id = %s
           AND g.biotype = 'protein_coding'
         GROUP BY g.gene_id, g.stable_id, t.stable_id, sr.name,
                  g.seq_region_start, g.seq_region_end
         ORDER BY COUNT(et.exon_id) DESC, g.gene_id
         LIMIT 1
    """
    with connection.cursor() as cursor:
        cursor.execute(query, (seq_region_id,))
        row = cursor.fetchone()
    return SampleGene(*row) if row else None


def choose_sample_gene(connection: pymysql.connections.Connection) -> SampleGene:
    """Apply the HiveAddPlaceholderLocation selection rules."""
    regions = longest_seq_regions(connection)
    if not regions:
        raise RuntimeError("No sequence regions were found in the core database")

    with connection.cursor() as cursor:
        cursor.execute("SELECT COUNT(*) FROM supporting_feature")
        has_supporting = cursor.fetchone()[0] > 0

    if has_supporting:
        for region_id in regions:
            sample_gene = select_supported_transcript(connection, region_id)
            if sample_gene:
                return sample_gene
    else:
        for region_id in regions:
            sample_gene = select_largest_gene(connection, region_id)
            if sample_gene:
                return sample_gene

    raise RuntimeError("Could not obtain a suitable sample gene and transcript")


def store_sample_metadata(
    connection: pymysql.connections.Connection, sample_gene: SampleGene, species_id: int
) -> None:
    """Replace the six sample metadata keys in a single transaction."""
    metadata = (
        ("sample.location_param", sample_gene.location),
        ("sample.location_text", sample_gene.location),
        ("sample.gene_param", sample_gene.gene_stable_id),
        ("sample.gene_text", sample_gene.gene_stable_id),
        ("sample.transcript_param", sample_gene.transcript_stable_id),
        ("sample.transcript_text", sample_gene.transcript_stable_id),
    )
    keys = [key for key, _ in metadata]
    placeholders = ", ".join(["%s"] * len(keys))

    with connection.cursor() as cursor:
        cursor.execute(
            f"DELETE FROM meta WHERE species_id = %s AND meta_key IN ({placeholders})",
            (species_id, *keys),
        )
        cursor.executemany(
            "INSERT INTO meta (species_id, meta_key, meta_value) VALUES (%s, %s, %s)",
            [(species_id, key, value) for key, value in metadata],
        )
    connection.commit()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-name", required=True, help="Core database name")
    parser.add_argument("--host", required=True, help="MySQL host")
    parser.add_argument("--port", required=True, type=int, help="MySQL port")
    parser.add_argument("--user", required=True, help="MySQL user")
    parser.add_argument("--password", default="", help="MySQL password")
    parser.add_argument("--species-id", type=int, default=1, help="Core species_id")
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    connection = pymysql.connect(
        host=args.host,
        port=args.port,
        user=args.user,
        password=args.password,
        database=args.db_name,
        autocommit=False,
    )
    try:
        sample_gene = choose_sample_gene(connection)
        store_sample_metadata(connection, sample_gene, args.species_id)
        LOGGER.info("Stored sample metadata for %s (%s)", sample_gene.gene_stable_id, sample_gene.location)
    finally:
        connection.close()


if __name__ == "__main__":
    main()
