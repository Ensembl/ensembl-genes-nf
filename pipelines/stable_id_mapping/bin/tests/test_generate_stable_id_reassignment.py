from __future__ import annotations

import json
import sys
from argparse import Namespace
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "bin"),
)

import generate_stable_id_reassignment
from stable_id_mapping.ids import parse_id_range


def make_args(
    tmp_path: Path,
    range_end: int = 199,
) -> Namespace:
    return Namespace(
        db_name="test_core",
        gene_range=parse_id_range(
            f"ENSXG:00000000100-{range_end:011d}"
        ),
        transcript_range=parse_id_range(
            f"ENSXT:00000000100-{range_end:011d}"
        ),
        translation_range=parse_id_range(
            f"ENSXP:00000000100-{range_end:011d}"
        ),
        output_sql=tmp_path / "reassign.sql",
        dry_run_sql=tmp_path / "dry_run.sql",
        output_json=tmp_path / "summary.json",
        batch_size=500,
        backup_prefix="test_backup",
    )


def empty_populations() -> dict[str, list[dict]]:
    return {
        "gene": [],
        "transcript": [],
        "translation": [],
        "exon": [],
    }


def test_all_correct_population_produces_noop(
    tmp_path: Path,
    monkeypatch,
) -> None:
    args = make_args(tmp_path)
    populations = {
        "gene": [
            {
                "feature_id": 1,
                "stable_id": "ENSXG00000000100",
                "version": 1,
            }
        ],
        "transcript": [
            {
                "feature_id": 2,
                "stable_id": "ENSXT00000000100",
                "version": 1,
            }
        ],
        "translation": [
            {
                "feature_id": 3,
                "stable_id": "ENSXP00000000100",
                "version": 1,
            }
        ],
        "exon": [
            {
                "feature_id": 4,
                "stable_id": "ENSXE00000000100",
                "version": 1,
            }
        ],
    }

    monkeypatch.setattr(
        generate_stable_id_reassignment,
        "parse_args",
        lambda: args,
    )
    monkeypatch.setattr(
        generate_stable_id_reassignment,
        "load_current_features",
        lambda db_name: populations,
    )

    assert generate_stable_id_reassignment.main() == 0

    summary = json.loads(args.output_json.read_text())
    assert summary["status"] == "clean"
    assert "No stable-ID reassignment is required" in (
        args.output_sql.read_text()
    )
    assert "UPDATE" not in args.output_sql.read_text()


def test_mixed_population_warns_and_reassigns_everything(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    args = make_args(tmp_path)
    populations = empty_populations()
    populations["gene"] = [
        {
            "feature_id": 1,
            "stable_id": "ENSXG00000000100",
            "version": 5,
        },
        {
            "feature_id": 2,
            "stable_id": "ENSXG00000000999",
            "version": 7,
        },
    ]

    monkeypatch.setattr(
        generate_stable_id_reassignment,
        "parse_args",
        lambda: args,
    )
    monkeypatch.setattr(
        generate_stable_id_reassignment,
        "load_current_features",
        lambda db_name: populations,
    )

    assert generate_stable_id_reassignment.main() == 0

    captured = capsys.readouterr()
    assert "MIXED STABLE-ID POPULATION" in captured.err

    summary = json.loads(args.output_json.read_text())
    assert summary["status"] == "reassignment_generated"
    assert summary["assignments"]["gene"] == 2

    sql = args.output_sql.read_text()
    assert "ENSXG00000000100" in sql
    assert "ENSXG00000000101" in sql
    assert "f.version = 1" in sql


def test_exhausted_range_fails_with_clear_message(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    args = make_args(tmp_path, range_end=100)
    populations = empty_populations()
    populations["gene"] = [
        {
            "feature_id": 1,
            "stable_id": "WRONG1",
            "version": 1,
        },
        {
            "feature_id": 2,
            "stable_id": "WRONG2",
            "version": 1,
        },
    ]

    monkeypatch.setattr(
        generate_stable_id_reassignment,
        "parse_args",
        lambda: args,
    )
    monkeypatch.setattr(
        generate_stable_id_reassignment,
        "load_current_features",
        lambda db_name: populations,
    )

    assert generate_stable_id_reassignment.main() == 2

    captured = capsys.readouterr()
    assert "REGISTRY RANGE EXHAUSTED" in captured.err
    assert "feature type: gene" in captured.err
    assert "required=2" in captured.err
    assert "available=1" in captured.err
    assert "ENSXG:00000000100-00000000100" in captured.err

    summary = json.loads(args.output_json.read_text())
    assert summary["status"] == "fatal_range_exhausted"
