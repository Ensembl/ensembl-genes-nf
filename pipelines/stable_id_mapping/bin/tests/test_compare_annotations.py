from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import compare_annotations
from stable_id_mapping.models import Decision


def test_comparison_writes_only_final_mapped_relationships(
    monkeypatch,
    tmp_path: Path,
) -> None:
    features = {
        "gene": {},
        "transcript": {},
        "translation": {},
        "exon": {},
    }
    monkeypatch.setattr(compare_annotations, "parse_gff3", lambda *args, **kwargs: features)
    monkeypatch.setattr(compare_annotations, "parse_missing_gene_ids", lambda path: set())
    monkeypatch.setattr(
        compare_annotations,
        "load_lifton_score_evidence",
        lambda *args: SimpleNamespace(),
    )
    monkeypatch.setattr(
        compare_annotations,
        "load_mapping_rules",
        lambda path: SimpleNamespace(
            coordinate_overlap=SimpleNamespace(min_overlap=0.75)
        ),
    )
    monkeypatch.setattr(
        compare_annotations,
        "build_decisions",
        lambda *args, **kwargs: [
            Decision(
                feature_type="gene",
                action="mapped",
                current_stable_id="TARGETG1",
                current_version=1,
                old_stable_id="REFG1",
                old_version=1,
                new_stable_id="REFG1",
                new_version=2,
                mapping_session_id=0,
                score=0.9,
                reason="mapped gene matched a target gene by lifton structural evidence",
            ),
            Decision(
                feature_type="gene",
                action="missing",
                current_stable_id=None,
                current_version=0,
                old_stable_id="REFG2",
                old_version=1,
                new_stable_id=None,
                new_version=0,
                mapping_session_id=0,
                score=0.0,
                reason="reference gene was not projected by LiftOn",
            ),
            Decision(
                feature_type="gene",
                action="new",
                current_stable_id="TARGETG2",
                current_version=1,
                old_stable_id=None,
                old_version=0,
                new_stable_id="COMPAREG000000000001",
                new_version=1,
                mapping_session_id=0,
                score=0.0,
                reason="target gene was not claimed by any mapped gene",
            ),
        ],
    )

    output_tsv = tmp_path / "relationships.tsv"
    compare_annotations.run_comparison(
        ref_gff=tmp_path / "ref.gff3",
        target_gff=tmp_path / "target.gff3",
        mapped_gff=tmp_path / "mapped.gff3",
        missing_report=tmp_path / "missing.txt",
        transcript_pairs=tmp_path / "transcript_pairs.tsv",
        gene_pairs=tmp_path / "gene_pairs.tsv",
        rules_config=tmp_path / "rules.json",
        output_tsv=output_tsv,
        include_translations=True,
    )

    assert output_tsv.read_text(encoding="utf-8").splitlines() == [
        "feature_type\treference_stable_id\ttarget_stable_id\tmethod\tscore\treason",
        "gene\tREFG1\tTARGETG1\tstructural\t0.900000\tmapped gene matched a target gene by lifton structural evidence",
    ]
