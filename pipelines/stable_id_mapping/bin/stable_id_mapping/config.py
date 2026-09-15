"""Configuration objects for stable-ID event generation."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

from .ids import StableIdRange
from .scoring import ScoreEvidence


def default_backup_prefix() -> str:
    return f"stable_id_mapper_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}"


@dataclass(frozen=True)
class StableIdEventConfig:
    ref_gff: Path
    target_gff: Path
    mapped_gff: Path
    report: Path
    mapping_session_id: int
    gene_range: StableIdRange
    transcript_range: StableIdRange
    translation_range: StableIdRange
    output_sql: Path
    ref_fasta: Optional[Path] = None
    target_fasta: Optional[Path] = None
    output_tsv: Optional[Path] = None
    db_name: Optional[str] = None
    include_translations: bool = False
    dry_run: bool = False
    backup_prefix: str = field(default_factory=default_backup_prefix)
    batch_size: int = 500
    replace_events_for_session: bool = False
    min_overlap: float = 0.75
    score_evidence: ScoreEvidence = field(default_factory=ScoreEvidence)
    output_score_evidence_tsv: Optional[Path] = None

    def validate(self) -> None:
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")
        if self.min_overlap < 0:
            raise ValueError("min_overlap must be >= 0")
        if self.ref_fasta is None or self.target_fasta is None:
            raise ValueError(
                "ref_fasta and target_fasta are required "
                "for stable-ID version comparison"
            )
        for path in (
            self.ref_gff,
            self.target_gff,
            self.mapped_gff,
            self.report,
            self.ref_fasta,
            self.target_fasta,
        ):
            if not path.exists():
                raise FileNotFoundError(path)
