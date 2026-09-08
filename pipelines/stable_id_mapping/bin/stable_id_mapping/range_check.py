"""Compare existing stable IDs with registry allocations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .ids import StableIdRange

from collections import Counter
from collections.abc import Iterable

import os

import pymysql


STABLE_ID_NUMERIC_WIDTH = 11

FEATURE_TABLES = {
    "gene": ("gene", "gene_id"),
    "transcript": ("transcript", "transcript_id"),
    "translation": ("translation", "translation_id"),
    "exon": ("exon", "exon_id"),
}

class DuplicateStableIdError(ValueError):
    """Raised when a stable ID occurs more than once."""


@dataclass(frozen=True)
class StableIdRangeCheck:
    agrees: bool
    reason: str
    numeric_value: Optional[int] = None


@dataclass(frozen=True)
class StableIdPopulationCheck:
    total: int
    agreeing: int
    disagreeing: int
    reason_counts: dict[str, int]
    checks: tuple[StableIdRangeCheck, ...]


def connect_read_only(db_name: str):
    return pymysql.connect(
        host=os.environ["GBS1"],
        port=int(os.environ["GBP1"]),
        user="ensro",
        password="",
        database=db_name,
        cursorclass=pymysql.cursors.DictCursor,
    )


def load_current_features(db_name: str) -> dict[str, list[dict]]:
    populations: dict[str, list[dict]] = {}
    connection = connect_read_only(db_name)

    try:
        with connection.cursor() as cursor:
            for feature_type, (table_name, primary_key) in FEATURE_TABLES.items():
                cursor.execute(
                    f"SELECT {primary_key} AS feature_id, "
                    "stable_id, COALESCE(version, 0) AS version "
                    f"FROM {table_name} "
                    f"ORDER BY {primary_key}"
                )
                populations[feature_type] = list(cursor.fetchall())
    finally:
        connection.close()

    return populations


def check_stable_id(
    stable_id: Optional[str],
    expected_range: StableIdRange,
) -> StableIdRangeCheck:
    """Check one stable ID against its registry-derived allocation."""

    if stable_id is None or stable_id == "":
        return StableIdRangeCheck(
            agrees=False,
            reason="missing",
        )

    if not stable_id.startswith(expected_range.prefix):
        return StableIdRangeCheck(
            agrees=False,
            reason="wrong_prefix",
        )

    numeric_text = stable_id[len(expected_range.prefix):]

    if not numeric_text.isdigit():
        return StableIdRangeCheck(
            agrees=False,
            reason="non_numeric",
        )

    if len(numeric_text) != STABLE_ID_NUMERIC_WIDTH:
        return StableIdRangeCheck(
            agrees=False,
            reason="wrong_width",
        )

    numeric_value = int(numeric_text)

    if not expected_range.start <= numeric_value <= expected_range.end:
        return StableIdRangeCheck(
            agrees=False,
            reason="outside_range",
            numeric_value=numeric_value,
        )

    return StableIdRangeCheck(
        agrees=True,
        reason="agrees",
        numeric_value=numeric_value,
    )

def find_duplicate_stable_ids(
    stable_ids: Iterable[Optional[str]],
) -> tuple[str, ...]:
    """Return every non-empty stable ID occurring more than once."""

    counts = Counter(
        stable_id
        for stable_id in stable_ids
        if stable_id is not None and stable_id != ""
    )

    return tuple(
        sorted(
            stable_id
            for stable_id, count in counts.items()
            if count > 1
        )
    )


def check_stable_id_population(
    stable_ids: Iterable[Optional[str]],
    expected_range: StableIdRange,
) -> StableIdPopulationCheck:
    """Compare a stable-ID population with one registry allocation."""

    stable_ids = tuple(stable_ids)
    duplicates = find_duplicate_stable_ids(stable_ids)

    if duplicates:
        raise DuplicateStableIdError(
            "Duplicate stable IDs detected: "
            + ", ".join(duplicates)
        )

    checks = tuple(
        check_stable_id(stable_id, expected_range)
        for stable_id in stable_ids
    )
    reason_counts = Counter(check.reason for check in checks)
    agreeing = reason_counts["agrees"]

    return StableIdPopulationCheck(
        total=len(checks),
        agreeing=agreeing,
        disagreeing=len(checks) - agreeing,
        reason_counts=dict(sorted(reason_counts.items())),
        checks=checks,
    )

