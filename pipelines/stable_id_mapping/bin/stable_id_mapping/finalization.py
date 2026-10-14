"""SQL generation for stable-ID-dependent sample-gene metadata."""

from __future__ import annotations

from typing import TextIO


SAMPLE_GENE_KEYS = (
    "sample.gene_param",
    "sample.gene_text",
    "genebuild.sample_gene",
)


def sql_string(value: str | None) -> str:
    if value is None:
        return "NULL"
    return "'" + value.replace("\\", "\\\\").replace("'", "''") + "'"


def write_metadata_map_table(
    handle: TextIO,
    rows: list[tuple[str, str, str]],
    *,
    batch_size: int,
) -> None:
    """Create the old/current-to-final stable-ID map used by sample metadata."""
    table = "tmp_stable_id_mapper_metadata_map"
    handle.write(
        f"CREATE TEMPORARY TABLE {table} (\n"
        "  feature_type VARCHAR(32) NOT NULL,\n"
        "  old_stable_id VARCHAR(128) NOT NULL,\n"
        "  new_stable_id VARCHAR(128) NOT NULL,\n"
        "  PRIMARY KEY (feature_type, old_stable_id)\n"
        ");\n\n"
    )
    values = [
        "(" + ", ".join(
            [sql_string(feature_type), sql_string(old_id), sql_string(new_id)]
        ) + ")"
        for feature_type, old_id, new_id in rows
    ]
    for start in range(0, len(values), batch_size):
        handle.write(
            f"INSERT INTO {table} (feature_type, old_stable_id, new_stable_id) VALUES\n"
            + ",\n".join(values[start : start + batch_size])
            + "\nON DUPLICATE KEY UPDATE new_stable_id = VALUES(new_stable_id);\n\n"
        )


def write_sample_gene_finalization_sql(handle: TextIO) -> None:
    """Update sample-gene metadata using the final gene stable IDs."""
    gene_keys = ", ".join(sql_string(key) for key in SAMPLE_GENE_KEYS)
    handle.write(
        "-- Remap stable-ID-dependent sample-gene metadata after gene updates.\n"
        "CREATE TEMPORARY TABLE tmp_stable_id_mapper_sample_gene_resolution AS\n"
        "SELECT m.species_id, m.meta_key,\n"
        "       COALESCE(map.new_stable_id, (\n"
        "         SELECT g.stable_id\n"
        "           FROM gene g\n"
        "           JOIN seq_region sr ON sr.seq_region_id = g.seq_region_id\n"
        "          WHERE sr.name = SUBSTRING_INDEX(loc.location_value, ':', 1)\n"
        "            AND g.seq_region_start <= CAST(SUBSTRING_INDEX(SUBSTRING_INDEX(loc.location_value, ':', -1), '-', -1) AS UNSIGNED)\n"
        "            AND g.seq_region_end >= CAST(SUBSTRING_INDEX(SUBSTRING_INDEX(loc.location_value, ':', -1), '-', 1) AS UNSIGNED)\n"
        "          ORDER BY ABS((g.seq_region_start + g.seq_region_end) -\n"
        "                       (CAST(SUBSTRING_INDEX(SUBSTRING_INDEX(loc.location_value, ':', -1), '-', 1) AS UNSIGNED) +\n"
        "                        CAST(SUBSTRING_INDEX(SUBSTRING_INDEX(loc.location_value, ':', -1), '-', -1) AS UNSIGNED))),\n"
        "                   g.gene_id\n"
        "          LIMIT 1\n"
        "       )) AS resolved_value\n"
        "  FROM meta m\n"
        "  LEFT JOIN tmp_stable_id_mapper_metadata_map map\n"
        "    ON map.feature_type = 'gene' AND map.old_stable_id = m.meta_value\n"
        "  LEFT JOIN (\n"
        "        SELECT species_id,\n"
        "               COALESCE(\n"
        "                 MAX(CASE WHEN meta_key = 'sample.location_param' THEN meta_value END),\n"
        "                 MAX(CASE WHEN meta_key = 'sample.location_text' THEN meta_value END)\n"
        "               ) AS location_value\n"
        "          FROM meta\n"
        "         WHERE meta_key IN ('sample.location_param', 'sample.location_text')\n"
        "         GROUP BY species_id\n"
        "  ) loc ON loc.species_id = m.species_id\n"
        f" WHERE m.meta_key IN ({gene_keys});\n\n"
        "UPDATE meta m\n"
        "JOIN tmp_stable_id_mapper_sample_gene_resolution r\n"
        "  ON r.species_id = m.species_id AND r.meta_key = m.meta_key\n"
        "SET m.meta_value = r.resolved_value\n"
        "WHERE r.resolved_value IS NOT NULL;\n\n"
        "DROP PROCEDURE IF EXISTS stable_id_mapper_sample_gene_validation;\n"
        "DELIMITER //\n"
        "CREATE PROCEDURE stable_id_mapper_sample_gene_validation()\n"
        "BEGIN\n"
        "  IF EXISTS (\n"
        "    SELECT 1 FROM meta m\n"
        f"     WHERE m.meta_key IN ({gene_keys})\n"
        "       AND NOT EXISTS (SELECT 1 FROM gene g WHERE g.stable_id = m.meta_value)\n"
        "  ) THEN\n"
        "    SIGNAL SQLSTATE '45000' SET MESSAGE_TEXT = 'stale sample gene stable ID in meta';\n"
        "  END IF;\n"
        "END//\n"
        "DELIMITER ;\n"
        "CALL stable_id_mapper_sample_gene_validation();\n"
        "DROP PROCEDURE stable_id_mapper_sample_gene_validation;\n\n"
    )
