#!/usr/bin/env python3
"""Shared I/O, normalized-record, filtering, and reporting helpers."""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple
from urllib.parse import unquote


NORMALIZED_FIELDS = (
    "seqid", "source", "feature_type", "start", "end", "score", "strand",
    "phase", "feature_id", "gene_id", "parent_ids", "transcript_ids",
    "biotype", "is_canonical",
)

CANONICAL_TAGS = {"canonical", "ensembl_canonical", "mane_select"}
BIOTYPE_KEYS = (
    "biotype", "gene_biotype", "transcript_biotype", "gene_type", "transcript_type",
)
TRUE_VALUES = {"1", "true", "yes"}
GTF_ATTRIBUTE = re.compile(r'(\S+)\s+"([^"]*)"')


@dataclass(frozen=True)
class Feature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    feature_id: str
    gene_id: str
    parent_ids: Tuple[str, ...]
    transcript_ids: Tuple[str, ...]
    biotype: str
    is_canonical: bool


def open_text(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def split_annotation_line(line: str) -> List[str]:
    stripped = line.rstrip("\r\n")
    if "\t" in stripped:
        return stripped.split("\t")
    return stripped.split(None, 8)


def parse_gff3_attributes(text: str) -> Dict[str, object]:
    result: Dict[str, object] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, value = item.partition("=")
        result[key.strip()] = unquote(value.strip())
    return result


def parse_gtf_attributes(text: str) -> Dict[str, object]:
    result: Dict[str, object] = {}
    for key, value in GTF_ATTRIBUTE.findall(text):
        if key == "tag":
            result.setdefault(key, [])
            assert isinstance(result[key], list)
            result[key].append(value)
        else:
            result[key] = value
    if result:
        return result
    for item in text.strip().strip(";").split(";"):
        parts = item.strip().split(None, 1)
        if len(parts) == 2:
            result[parts[0]] = parts[1].strip().strip('"')
    return result


def first(attrs: Dict[str, object], keys: Sequence[str]) -> str:
    for key in keys:
        value = attrs.get(key, "")
        if isinstance(value, list):
            if value:
                return str(value[0])
        elif value:
            return str(value)
    return ""


def values(attrs: Dict[str, object], keys: Sequence[str]) -> Tuple[str, ...]:
    value = first(attrs, keys)
    return tuple(part.strip().strip('"') for part in value.split(",") if part.strip())


def is_canonical(attrs: Dict[str, object]) -> bool:
    tags = attrs.get("tag", [])
    if isinstance(tags, str):
        tags = [tags]
    if any(str(tag).strip().lower() in CANONICAL_TAGS for tag in tags):
        return True
    return any(
        str(attrs.get(key, "")).strip().strip('"').lower() in TRUE_VALUES
        for key in ("is_canonical", "transcript_is_canonical")
    )


def parse_biotypes(values_arg: Optional[Sequence[str]]) -> Optional[Set[str]]:
    if not values_arg:
        return None
    result = set()
    for value in values_arg:
        result.update(item.strip().lower() for item in value.split(",") if item.strip())
    return result or None


def feature_matches_biotype(feature: Feature, allowed: Optional[Set[str]]) -> bool:
    return not allowed or feature.biotype.lower() in allowed


def resolve_selected_features(
    features: Sequence[Feature],
    include_biotypes: Optional[Set[str]] = None,
    exclude_biotypes: Optional[Set[str]] = None,
    canonical_only: bool = False,
) -> List[Feature]:
    """Apply filters after A has read the annotation once.

    Canonical status is inherited from directly canonical transcript/gene
    records, preserving the two-pass semantic of the original scripts while
    rereading only the compact normalized records.
    """
    canonical_transcripts = {
        feature.feature_id
        for feature in features
        if feature.is_canonical and feature.feature_type.lower() in {
            "transcript", "mrna", "rrna", "trna", "ncrna", "lncrna", "snrna", "snorna"
        }
    }
    canonical_genes = {
        feature.feature_id
        for feature in features
        if feature.is_canonical and feature.feature_type.lower() == "gene"
    }
    canonical_transcripts.update(
        feature.feature_id
        for feature in features
        if feature.feature_id
        and set(feature.parent_ids) & canonical_genes
        and feature.feature_type.lower() in {
            "transcript", "mrna", "rrna", "trna", "ncrna", "lncrna", "snrna", "snorna"
        }
    )
    canonical_genes.update(
        parent_id
        for feature in features
        if feature.feature_id in canonical_transcripts
        for parent_id in feature.parent_ids
    )

    included_ids = set()
    excluded_ids = set()
    if include_biotypes:
        for feature in features:
            if feature.biotype.lower() in include_biotypes:
                included_ids.update((feature.feature_id, feature.gene_id))
                included_ids.update(feature.parent_ids)
                included_ids.update(feature.transcript_ids)
    if exclude_biotypes:
        for feature in features:
            if feature.biotype.lower() in exclude_biotypes:
                excluded_ids.update((feature.feature_id, feature.gene_id))
                excluded_ids.update(feature.parent_ids)
                excluded_ids.update(feature.transcript_ids)

    selected = []
    for feature in features:
        related_ids = {feature.feature_id, feature.gene_id} | set(feature.parent_ids) | set(feature.transcript_ids)
        if include_biotypes and feature.biotype.lower() not in include_biotypes and not related_ids & included_ids:
            continue
        if exclude_biotypes and (feature.biotype.lower() in exclude_biotypes or related_ids & excluded_ids):
            continue
        if canonical_only:
            transcript_match = bool(set(feature.transcript_ids) & canonical_transcripts)
            parent_match = bool(set(feature.parent_ids) & canonical_transcripts)
            gene_match = feature.gene_id in canonical_genes or feature.feature_id in canonical_genes
            direct_match = feature.is_canonical
            if not (transcript_match or parent_match or gene_match or direct_match):
                continue
        selected.append(feature)
    return selected


def write_metrics(rows: Sequence[Tuple[str, object]], path: str) -> None:
    with open_text(path, "wt") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("metrics_name", "metrics_value"))
        for name, value in rows:
            writer.writerow((name, value))


def read_normalized(path: str) -> List[Feature]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(NORMALIZED_FIELDS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Normalized file is missing columns: " + ", ".join(sorted(missing)))
        result = []
        for row in reader:
            result.append(
                Feature(
                    seqid=row["seqid"], source=row["source"], feature_type=row["feature_type"],
                    start=int(row["start"]), end=int(row["end"]), score=row["score"],
                    strand=row["strand"], phase=row["phase"], feature_id=row["feature_id"],
                    gene_id=row["gene_id"],
                    parent_ids=tuple(x for x in row["parent_ids"].split(",") if x),
                    transcript_ids=tuple(x for x in row["transcript_ids"].split(",") if x),
                    biotype=row["biotype"], is_canonical=row["is_canonical"] == "1",
                )
            )
    return result

