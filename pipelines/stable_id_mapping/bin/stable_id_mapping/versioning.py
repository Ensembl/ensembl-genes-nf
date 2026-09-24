"""Feature-specific stable-ID version rules."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

import pysam
from Bio.Seq import Seq

from .gff3 import TRANSCRIPT_FEATURE_TYPES, parent_ids, parse_attrs, split_stable_id
from .models import Decision, Feature


@dataclass(frozen=True)
class CdsBlock:
    start: int
    end: int
    phase: int


@dataclass(frozen=True)
class AnnotationDetails:
    transcript_exons: dict[str, tuple[tuple[int, int], ...]]
    transcript_cds: dict[str, tuple[tuple[int, int], ...]]
    translation_cds: dict[str, tuple[CdsBlock, ...]]
    translation_tables: dict[str, int]


def parse_annotation_details(path: Path) -> AnnotationDetails:
    transcript_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    transcript_cds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    translation_cds: dict[str, list[CdsBlock]] = defaultdict(list)
    translation_tables: dict[str, int] = {}

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            (
                _seqid,
                _source,
                feature_type,
                start,
                end,
                _score,
                _strand,
                phase_text,
                attrs_text,
            ) = fields

            feature_type_lc = feature_type.lower()
            attrs = parse_attrs(attrs_text)
            parents = parent_ids(attrs.get("Parent"))
            start_i = int(start)
            end_i = int(end)

            if feature_type_lc == "exon":
                for parent_stable_id, _parent_version in parents:
                    transcript_exons[parent_stable_id].append(
                        (start_i, end_i)
                    )
                continue

            if feature_type_lc != "cds":
                continue

            phase = (
                int(phase_text)
                if phase_text in {"0", "1", "2"}
                else 0
            )

            protein_stable_id, _protein_version = split_stable_id(
                attrs.get("protein_id")
            )

            translation_table_text = (
                attrs.get("transl_table")
                or attrs.get("translation_table")
                or "1"
            )

            try:
                translation_table = int(translation_table_text)
            except ValueError:
                translation_table = 1

            for parent_stable_id, _parent_version in parents:
                transcript_cds[parent_stable_id].append(
                    (start_i, end_i)
                )

                translation_stable_id = (
                    protein_stable_id or parent_stable_id
                )

                translation_cds[translation_stable_id].append(
                    CdsBlock(
                        start=start_i,
                        end=end_i,
                        phase=phase,
                    )
                )

                previous_table = translation_tables.get(
                    translation_stable_id
                )
                if (
                    previous_table is not None
                    and previous_table != translation_table
                ):
                    raise ValueError(
                        f"{path}: translation {translation_stable_id} "
                        "uses conflicting translation tables"
                    )

                translation_tables[translation_stable_id] = (
                    translation_table
                )

    return AnnotationDetails(
        transcript_exons={
            stable_id: tuple(sorted(blocks))
            for stable_id, blocks in transcript_exons.items()
        },
        transcript_cds={
            stable_id: tuple(sorted(blocks))
            for stable_id, blocks in transcript_cds.items()
        },
        translation_cds={
            stable_id: tuple(
                sorted(
                    blocks,
                    key=lambda block: (block.start, block.end),
                )
            )
            for stable_id, blocks in translation_cds.items()
        },
        translation_tables=translation_tables,
    )


def ensure_fasta_index(path: Path) -> None:
    index_path = Path(f"{path}.fai")
    if not index_path.exists():
        pysam.faidx(str(path))


def ordered_blocks(
    blocks: tuple[tuple[int, int], ...],
    strand: str,
) -> tuple[tuple[int, int], ...]:
    return tuple(
        sorted(
            blocks,
            key=lambda block: block[0],
            reverse=strand == "-",
        )
    )


def sequence_from_blocks(
    fasta: pysam.FastaFile,
    feature: Feature,
    blocks: tuple[tuple[int, int], ...],
) -> str:
    parts: list[str] = []

    for start, end in ordered_blocks(blocks, feature.strand):
        sequence = fasta.fetch(
            feature.seqid,
            start - 1,
            end,
        )

        if feature.strand == "-":
            sequence = str(Seq(sequence).reverse_complement())

        parts.append(sequence.upper())

    return "".join(parts)


def transcript_blocks(
    details: AnnotationDetails,
    stable_id: str,
) -> tuple[tuple[int, int], ...]:
    return (
        details.transcript_exons.get(stable_id)
        or details.transcript_cds.get(stable_id)
        or ()
    )


def normalized_splice_pattern(
    feature: Feature,
    details: AnnotationDetails,
) -> tuple[tuple[int, int], ...]:
    blocks = transcript_blocks(details, feature.stable_id)
    if not blocks:
        return ()

    ordered = ordered_blocks(blocks, feature.strand)

    if feature.strand == "-":
        anchor = ordered[0][1]
        return tuple(
            (anchor - end, anchor - start)
            for start, end in ordered
        )

    anchor = ordered[0][0]
    return tuple(
        (start - anchor, end - anchor)
        for start, end in ordered
    )


def transcript_cdna(
    fasta: pysam.FastaFile,
    feature: Feature,
    details: AnnotationDetails,
) -> Optional[str]:
    blocks = details.transcript_exons.get(feature.stable_id)
    if not blocks:
        return None

    return sequence_from_blocks(
        fasta,
        feature,
        blocks,
    )


def translation_peptide(
    fasta: pysam.FastaFile,
    feature: Feature,
    details: AnnotationDetails,
) -> Optional[str]:
    cds_blocks = details.translation_cds.get(feature.stable_id)
    if not cds_blocks:
        return None

    ordered = tuple(
        sorted(
            cds_blocks,
            key=lambda block: block.start,
            reverse=feature.strand == "-",
        )
    )

    coordinate_blocks = tuple(
        (block.start, block.end)
        for block in ordered
    )

    cds_sequence = sequence_from_blocks(
        fasta,
        feature,
        coordinate_blocks,
    )

    initial_phase = ordered[0].phase
    if initial_phase:
        cds_sequence = cds_sequence[initial_phase:]

    complete_length = len(cds_sequence) - (len(cds_sequence) % 3)
    cds_sequence = cds_sequence[:complete_length]

    if not cds_sequence:
        return ""

    translation_table = details.translation_tables.get(
        feature.stable_id,
        1,
    )

    peptide = str(
        Seq(cds_sequence).translate(
            table=translation_table,
            to_stop=False,
        )
    )

    return peptide.removesuffix("*")


def transcript_membership_by_gene(
    transcripts: dict[str, Feature],
) -> dict[str, set[str]]:
    memberships: dict[str, set[str]] = defaultdict(set)

    for transcript in transcripts.values():
        if transcript.parent_stable_id:
            memberships[transcript.parent_stable_id].add(
                transcript.stable_id
            )

    return memberships


def apply_version_rules(
    decisions: list[Decision],
    ref_features: dict[str, dict[str, Feature]],
    target_features: dict[str, dict[str, Feature]],
    mapped_features: dict[str, dict[str, Feature]],
    ref_gff: Path,
    target_gff: Path,
    mapped_gff: Path,
    ref_fasta: Path,
    target_fasta: Path,
) -> list[Decision]:
    ref_details = parse_annotation_details(ref_gff)
    target_details = parse_annotation_details(target_gff)
    mapped_details = parse_annotation_details(mapped_gff)

    ref_memberships = transcript_membership_by_gene(
        ref_features["transcript"]
    )

    assigned_transcript_ids = {
        decision.current_stable_id: decision.new_stable_id
        for decision in decisions
        if decision.feature_type == "transcript"
        and decision.current_stable_id
        and decision.new_stable_id
    }

    target_memberships: dict[str, set[str]] = defaultdict(set)
    for transcript in target_features["transcript"].values():
        if not transcript.parent_stable_id:
            continue

        assigned_id = assigned_transcript_ids.get(
            transcript.stable_id
        )
        if assigned_id:
            target_memberships[
                transcript.parent_stable_id
            ].add(assigned_id)

    ensure_fasta_index(ref_fasta)
    ensure_fasta_index(target_fasta)

    updated: list[Decision] = []

    with pysam.FastaFile(str(ref_fasta)) as ref_sequences:
        with pysam.FastaFile(str(target_fasta)) as target_sequences:
            for decision in decisions:
                if decision.action != "mapped":
                    updated.append(decision)
                    continue

                changes: list[str] = []

                if decision.feature_type == "gene":
                    old_transcripts = ref_memberships.get(
                        decision.old_stable_id or "",
                        set(),
                    )
                    new_transcripts = target_memberships.get(
                        decision.current_stable_id or "",
                        set(),
                    )

                    if old_transcripts != new_transcripts:
                        changes.append("transcript set changed")

                elif decision.feature_type == "transcript":
                    old = ref_features["transcript"].get(
                        decision.old_stable_id or ""
                    )
                    target = target_features["transcript"].get(
                        decision.current_stable_id or ""
                    )
                    projected = mapped_features["transcript"].get(
                        decision.old_stable_id or ""
                    )

                    if old is None or target is None:
                        changes.append("transcript comparison unavailable")
                    else:
                        if projected is None:
                            changes.append(
                                "projected transcript unavailable"
                            )
                        else:
                            old_pattern = normalized_splice_pattern(
                                projected,
                                mapped_details,
                            )
                            target_pattern = normalized_splice_pattern(
                                target,
                                target_details,
                            )

                            if old_pattern != target_pattern:
                                changes.append(
                                    "splicing pattern changed"
                                )

                            old_location = (
                                projected.seqid,
                                projected.start,
                                projected.end,
                                projected.strand,
                            )
                            target_location = (
                                target.seqid,
                                target.start,
                                target.end,
                                target.strand,
                            )

                            if old_location != target_location:
                                changes.append(
                                    "chromosome location changed"
                                )

                        old_cdna = transcript_cdna(
                            ref_sequences,
                            old,
                            ref_details,
                        )
                        target_cdna = transcript_cdna(
                            target_sequences,
                            target,
                            target_details,
                        )

                        if old_cdna is None or target_cdna is None:
                            changes.append(
                                "cDNA comparison unavailable"
                            )
                        elif old_cdna != target_cdna:
                            changes.append("cDNA sequence changed")

                elif decision.feature_type == "translation":
                    old = ref_features["translation"].get(
                        decision.old_stable_id or ""
                    )
                    target = target_features["translation"].get(
                        decision.current_stable_id or ""
                    )

                    if old is None or target is None:
                        changes.append("peptide comparison unavailable")
                    else:
                        old_peptide = translation_peptide(
                            ref_sequences,
                            old,
                            ref_details,
                        )
                        target_peptide = translation_peptide(
                            target_sequences,
                            target,
                            target_details,
                        )

                        if (
                            old_peptide is None
                            or target_peptide is None
                        ):
                            changes.append(
                                "peptide comparison unavailable"
                            )
                        elif old_peptide != target_peptide:
                            changes.append("peptide sequence changed")

                if changes:
                    new_version = decision.old_version + 1
                    version_reason = (
                        "version incremented: " + ", ".join(changes)
                    )
                else:
                    new_version = decision.old_version
                    version_reason = (
                        "version unchanged: relevant features unchanged"
                    )

                updated.append(
                    replace(
                        decision,
                        new_version=new_version,
                        reason=(
                            f"{decision.reason}; {version_reason}"
                        ),
                    )
                )

    return updated
