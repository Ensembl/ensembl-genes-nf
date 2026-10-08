from __future__ import annotations

from io import StringIO

from stable_id_mapping.finalization import (
    write_metadata_map_table,
    write_sample_gene_finalization_sql,
)


def test_metadata_map_contains_gene_and_transcript_aliases() -> None:
    handle = StringIO()

    write_metadata_map_table(
        handle,
        [
            ("gene", "OLDG", "NEWG"),
            ("transcript", "OLDT", "NEWT"),
        ],
        batch_size=500,
    )

    sql = handle.getvalue()
    assert "tmp_stable_id_mapper_metadata_map" in sql
    assert "'OLDG', 'NEWG'" in sql
    assert "'OLDT', 'NEWT'" in sql
    assert "ON DUPLICATE KEY UPDATE" in sql


def test_sample_gene_finalization_uses_location_fallback_and_hard_gate() -> None:
    handle = StringIO()

    write_sample_gene_finalization_sql(handle)

    sql = handle.getvalue()
    for key in (
        "sample.gene_param",
        "sample.gene_text",
        "genebuild.sample_gene",
    ):
        assert key in sql
    assert "sample.location_param" in sql
    assert "ORDER BY ABS" in sql
    assert "stale sample gene stable ID in meta" in sql
