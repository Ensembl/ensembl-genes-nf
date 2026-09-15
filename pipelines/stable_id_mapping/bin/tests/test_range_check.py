from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "bin"),
)

from stable_id_mapping.ids import (
    exon_range_from_gene_range,
    parse_id_range,
)
from stable_id_mapping.range_check import (
    check_stable_id,
    DuplicateStableIdError,
    check_stable_id_population,
    find_duplicate_stable_ids,
)


@pytest.mark.parametrize(
    ("range_text", "stable_id"),
    [
        ("ENSXG:00000000100-00000000199", "ENSXG00000000100"),
        ("ENSXT:00000000100-00000000199", "ENSXT00000000150"),
        ("ENSXP:00000000100-00000000199", "ENSXP00000000199"),
    ],
)
def test_stable_id_agrees_with_registry_range(
    range_text: str,
    stable_id: str,
) -> None:
    result = check_stable_id(
        stable_id,
        parse_id_range(range_text),
    )

    assert result.agrees is True
    assert result.reason == "agrees"


def test_exons_follow_the_same_registry_rules() -> None:
    gene_range = parse_id_range(
        "ENSXG:00000000100-00000000199"
    )
    exon_range = exon_range_from_gene_range(gene_range)

    result = check_stable_id(
        "ENSXE00000000125",
        exon_range,
    )

    assert result.agrees is True
    assert result.reason == "agrees"


@pytest.mark.parametrize(
    ("stable_id", "expected_reason"),
    [
        (None, "missing"),
        ("", "missing"),
        ("WRONG00000000150", "wrong_prefix"),
        ("ENSXGABCDEFGHIJK", "non_numeric"),
        ("ENSXG0000000150", "wrong_width"),
        ("ENSXG00000000099", "outside_range"),
        ("ENSXG00000000200", "outside_range"),
    ],
)
def test_stable_id_disagrees_with_registry_range(
    stable_id: str | None,
    expected_reason: str,
) -> None:
    registry_range = parse_id_range(
        "ENSXG:00000000100-00000000199"
    )

    result = check_stable_id(
        stable_id,
        registry_range,
    )

    assert result.agrees is False
    assert result.reason == expected_reason

def test_population_reports_agreement_and_disagreement() -> None:
    registry_range = parse_id_range(
        "ENSXG:00000000100-00000000199"
    )

    result = check_stable_id_population(
        [
            "ENSXG00000000100",
            "ENSXG00000000150",
            "ENSXG00000000200",
            "WRONG00000000125",
            None,
        ],
        registry_range,
    )

    assert result.total == 5
    assert result.agreeing == 2
    assert result.disagreeing == 3
    assert result.reason_counts == {
        "agrees": 2,
        "missing": 1,
        "outside_range": 1,
        "wrong_prefix": 1,
    }


def test_duplicate_stable_ids_are_fatal() -> None:
    registry_range = parse_id_range(
        "ENSXG:00000000100-00000000199"
    )

    with pytest.raises(
        DuplicateStableIdError,
        match="ENSXG00000000125",
    ):
        check_stable_id_population(
            [
                "ENSXG00000000125",
                "ENSXG00000000125",
            ],
            registry_range,
        )


def test_duplicate_check_can_cover_all_feature_types() -> None:
    duplicates = find_duplicate_stable_ids(
        [
            "ENSXG00000000125",
            "ENSXT00000000126",
            "ENSXP00000000127",
            "ENSXE00000000128",
            "ENSXG00000000125",
        ]
    )

    assert duplicates == ("ENSXG00000000125",)


def test_multiple_missing_ids_are_not_duplicates() -> None:
    assert find_duplicate_stable_ids(
        [None, None, "", ""]
    ) == ()
