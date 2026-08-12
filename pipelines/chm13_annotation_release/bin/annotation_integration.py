#!/usr/bin/env python3
"""Compare and integrate projected and manually curated GFF3 annotations.

The comparison is deliberately transcript-centric.  A decision table can
override the default policy for individual manual transcripts:

    manual_transcript_id,action
    ENST...,add_manual_preserved
    ENST...,keep_manual
    ENST...,keep_projected
    ENST...,skip
    ENST...,promote_canonical
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
TRANSCRIPT_TYPES = {
    "transcript", "mRNA", "lnc_RNA", "pseudogenic_transcript",
    "unconfirmed_transcript", "ncRNA", "miRNA", "snRNA", "snoRNA",
    "scRNA", "rRNA",
}
VALID_ACTIONS = {"add_manual_preserved", "keep_both", "keep_manual", "keep_projected", "skip", "promote_canonical", "mark_manual_curated"}


@dataclass
class Record:
    line: str
    feature: str
    attrs: dict[str, str]
    seqname: str
    start: int
    end: int
    strand: str


@dataclass
class Model:
    ident: str
    gene_id: str
    seqname: str
    start: int
    end: int
    strand: str
    feature: str
    records: list[Record] = field(default_factory=list)
    exons: list[tuple[int, int]] = field(default_factory=list)
    tags: set[str] = field(default_factory=set)
    cds: list[tuple[int, int]] = field(default_factory=list)
    five_utr: list[tuple[int, int]] = field(default_factory=list)
    three_utr: list[tuple[int, int]] = field(default_factory=list)


def parse_attrs(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            result[key] = val
    return result


def read_gff(path: Path, retain_records: bool = True) -> tuple[list[str], dict[str, list[Record]], dict[str, Model]]:
    headers: list[str] = []
    genes: dict[str, list[Record]] = defaultdict(list)
    transcripts: dict[str, Model] = {}
    with path.open() as handle:
        for raw in handle:
            if raw.startswith("#"):
                headers.append(raw)
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_attrs(fields[8])
            record = Record(raw, fields[2], attrs, fields[0], int(fields[3]), int(fields[4]), fields[6])
            if record.feature in GENE_TYPES:
                ident = attrs.get("ID", "").removeprefix("gene:")
                genes[ident].append(record)
            elif record.feature in TRANSCRIPT_TYPES:
                ident = attrs.get("ID", "").removeprefix("transcript:")
                gene_id = attrs.get("Parent", "").removeprefix("gene:")
                transcripts[ident] = Model(
                    ident, gene_id, record.seqname, record.start, record.end,
                    record.strand, record.feature, [record] if retain_records else [], [],
                    set(attrs.get("tag", "").split(",")),
                )
            elif record.feature == "exon":
                ident = attrs.get("Parent", "").removeprefix("transcript:")
                if ident in transcripts:
                    if retain_records:
                        transcripts[ident].records.append(record)
                    transcripts[ident].exons.append((record.start, record.end))
            elif record.feature in {"CDS", "five_prime_UTR", "three_prime_UTR"}:
                ident = attrs.get("Parent", "").removeprefix("transcript:")
                if ident in transcripts:
                    if retain_records:
                        transcripts[ident].records.append(record)
                    interval = (record.start, record.end)
                    if record.feature == "CDS":
                        transcripts[ident].cds.append(interval)
                    elif record.feature == "five_prime_UTR":
                        transcripts[ident].five_utr.append(interval)
                    else:
                        transcripts[ident].three_utr.append(interval)
    for model in transcripts.values():
        model.exons.sort()
    return headers, genes, transcripts


def overlaps(a: Model, b: Model) -> bool:
    return a.seqname == b.seqname and a.strand == b.strand and a.start <= b.end and b.start <= a.end


def classify(manual: Model, projected_genes: Iterable[Model], projected_tx: Iterable[Model]) -> tuple[str, list[str], list[str]]:
    if "readthrough_transcript" in manual.tags:
        return "readthrough_excluded", [], []
    genes = [g for g in projected_genes if overlaps(manual, g)]
    gene_ids = sorted({g.gene_id for g in genes})
    candidates = [t for t in projected_tx if any(t.gene_id == gid for gid in gene_ids)]
    tx_ids = []
    scored: list[tuple[int, int, Model]] = []
    for candidate in candidates:
        if "readthrough_transcript" in candidate.tags:
            continue
        if candidate.strand != manual.strand:
            continue
        span_overlap = max(0, min(manual.end, candidate.end) - max(manual.start, candidate.start) + 1)
        if span_overlap:
            tx_ids.append(candidate.ident)
        score = 0
        if manual.exons and candidate.exons and manual.exons == candidate.exons:
            score = 3
        elif (manual.exons and candidate.exons and len(manual.exons) == len(candidate.exons)
              and manual.exons[1:-1] == candidate.exons[1:-1]):
            score = 2
        elif span_overlap / max(1, manual.end - manual.start + 1) >= 0.8:
            score = 1
        if span_overlap:
            scored.append((score, -abs(manual.start - candidate.start) - abs(manual.end - candidate.end), candidate))
    tx_ids = sorted(set(tx_ids))
    if not genes:
        return "novel_locus", gene_ids, tx_ids
    if not scored:
        return "no_transcript_overlap", gene_ids, tx_ids
    best_score, _, _ = max(scored, key=lambda item: (item[0], item[1]))
    if best_score == 3:
        return "exact_exon_structure", gene_ids, tx_ids
    if best_score == 2:
        return "same_internal_exon_chain_UTR_difference", gene_ids, tx_ids
    if best_score == 1:
        return "same_locus_different_exon_chain", gene_ids, tx_ids
    return "overlapping_different_exon_chain", gene_ids, tx_ids


def interval_span(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in intervals)


def utr_relation(manual: Model, projected: Model) -> str:
    def relation(manual_parts, projected_parts, label):
        m = interval_span(manual_parts)
        p = interval_span(projected_parts)
        if m == p and manual_parts == projected_parts:
            return f"{label}=identical"
        if not m and not p:
            return f"{label}=absent_both"
        if not m:
            return f"{label}=manual_absent"
        if not p:
            return f"{label}=manual_only"
        if m > p:
            return f"{label}=manual_longer_by_{m-p}bp"
        if m < p:
            return f"{label}=manual_shorter_by_{p-m}bp"
        return f"{label}=same_length_different_boundaries"
    return ";".join((
        relation(manual.five_utr, projected.five_utr, "five_prime_UTR"),
        relation(manual.three_utr, projected.three_utr, "three_prime_UTR"),
    ))


def model_explanation(manual: Model, projected: Model) -> str:
    if manual.exons == projected.exons:
        exon_relation = "exons=identical"
    elif len(manual.exons) == len(projected.exons) and manual.exons[1:-1] == projected.exons[1:-1]:
        exon_relation = "exons=same_internal_chain_UTR_or_terminal_difference"
    elif set(manual.exons).issubset(set(projected.exons)):
        exon_relation = "exons=manual_subset"
    elif set(projected.exons).issubset(set(manual.exons)):
        exon_relation = "exons=manual_superset"
    else:
        exon_relation = "exons=different_chain"
    return ";".join((exon_relation, utr_relation(manual, projected)))


def model_comparison(manual: Model, projected: Model) -> dict[str, str]:
    return {
        "best_projected_transcript_id": projected.ident,
        "manual_exon_count": str(len(manual.exons)),
        "projected_exon_count": str(len(projected.exons)),
        "manual_transcript_start": str(manual.start),
        "manual_transcript_end": str(manual.end),
        "projected_transcript_start": str(projected.start),
        "projected_transcript_end": str(projected.end),
        "manual_five_prime_utr_bp": str(interval_span(manual.five_utr)),
        "projected_five_prime_utr_bp": str(interval_span(projected.five_utr)),
        "manual_three_prime_utr_bp": str(interval_span(manual.three_utr)),
        "projected_three_prime_utr_bp": str(interval_span(projected.three_utr)),
    }


def decision_annotation(category: str, projected_gene_ids: list[str], projected_tx_ids: list[str]) -> tuple[str, str, str]:
    """Return (review_required, reason, recommended_action)."""
    if category == "novel_locus":
        return "false", "no_projected_locus", "add_as_new_manual_gene"
    if category == "exact_exon_structure":
        return "false", "exact_transcript_model_match", "add_manual_transcript_to_existing_gene"
    if category == "same_internal_exon_chain_UTR_difference":
        return "true", "same_CDS_or_internal_exon_chain_but_UTR_or_terminal_difference", "confirm_parent_gene_and_retain_manual_isoform"
    if category == "readthrough_excluded":
        return "true", "readthrough_transcript_excluded_from_automatic_assignment", "HAVANA_decision_required"
    if len(projected_gene_ids) > 1:
        return "true", "manual_model_overlaps_multiple_projected_genes", "confirm_parent_gene"
    if category == "same_locus_different_exon_chain":
        return "true", "same_locus_but_different_exon_chain", "confirm_same_gene_or_separate_gene"
    if category == "overlapping_different_exon_chain":
        return "true", "overlapping_locus_but_different_exon_chain", "confirm_same_gene_or_separate_gene"
    if category == "no_transcript_overlap":
        return "true", "projected_gene_overlap_without_matching_transcript", "confirm_parent_gene"
    return "true", category, "HAVANA_decision_required"


def build_overlap_index(models: Iterable[Model]) -> dict[tuple[str, str], list[Model]]:
    index: dict[tuple[str, str], list[Model]] = defaultdict(list)
    for model in models:
        index[(model.seqname, model.strand)].append(model)
    for bucket in index.values():
        bucket.sort(key=lambda model: model.start)
    return index


def overlapping_models(model: Model, index: dict[tuple[str, str], list[Model]]) -> list[Model]:
    # The projected annotation is sorted by coordinate. Stop at the first
    # model that begins after the manual model ends.
    return [candidate for candidate in index.get((model.seqname, model.strand), [])
            if candidate.start <= model.end and candidate.end >= model.start]


def load_decisions(path: Path | None) -> dict[str, str]:
    if not path:
        return {}
    with path.open(newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first else ","
        decisions = {row["manual_transcript_id"]: row["action"] for row in csv.DictReader(handle, delimiter=delimiter)}
    invalid = sorted({action for action in decisions.values() if action not in VALID_ACTIONS})
    if invalid:
        raise ValueError(f"Unsupported decision action(s): {', '.join(invalid)}")
    return decisions


def add_source(record: Record, promotion_status: str = "review", metadata: dict[str, str] | None = None) -> str:
    fields = record.line.rstrip("\n").split("\t")
    fields[1] = "HAVANA_manual"
    fields[8] += (
        ";integration_source=HAVANA_manual"
        ";manual_source=HAVANA"
        ";manual_annotation_priority=gold"
        ";manual_canonical_candidate=true"
        f";manual_promotion_status={promotion_status}"
    )
    for key, value in (metadata or {}).items():
        fields[8] += f";{key}={value}"
    if promotion_status == "promoted":
        fields[8] += ";manual_canonical=true"
    return "\t".join(fields) + "\n"


def add_manual_curated(record_line: str, manual_transcript_ids: list[str]) -> str:
    fields = record_line.rstrip("\n").split("\t")
    fields[8] += (
        ";manual_curated=true;manual_source=HAVANA"
        ";manual_annotation_priority=gold"
        ";manual_curation_status=confirmed_exact_model"
        f";manual_evidence_transcripts={','.join(manual_transcript_ids)}"
    )
    return "\t".join(fields) + "\n"


def write_merged(projected: Path | None, manual: Path | None, manual_genes: dict[str, list[Record]], manual_models: dict[str, Model], decisions: dict[str, str], drop_projected: set[str], curated_projected: dict[str, list[str]], manual_metadata: dict[str, dict[str, str]], gene_metadata: dict[str, dict[str, str]], output: Path) -> None:
    with output.open("w") as out:
        if projected and projected.exists():
            with projected.open() as handle:
                for raw in handle:
                    if raw.startswith("#"):
                        out.write(raw)
                        continue
                    attrs = parse_attrs(raw.rstrip("\n").split("\t", 8)[8])
                    ident = attrs.get("ID", "").removeprefix("transcript:")
                    parent = attrs.get("Parent", "").removeprefix("transcript:")
                    if ident in drop_projected or parent in drop_projected:
                        continue
                    evidence = curated_projected.get(ident) or curated_projected.get(parent)
                    out.write(add_manual_curated(raw, evidence)) if evidence else out.write(raw)
        if manual:
            kept_gene_ids = {
                model.gene_id for model in manual_models.values()
                if decisions.get(model.ident) not in {"keep_projected", "skip", "mark_manual_curated"}
            }
            for gene_id in sorted(kept_gene_ids):
                for record in manual_genes.get(gene_id, []):
                    out.write(add_source(record, metadata=gene_metadata.get(gene_id)))
            for model in manual_models.values():
                action = decisions.get(model.ident)
                if action in {"keep_projected", "skip", "mark_manual_curated"}:
                    continue
                promotion_status = "promoted" if action == "promote_canonical" else "review"
                for record in model.records:
                    out.write(add_source(record, promotion_status, manual_metadata.get(model.ident)))


def write_session(output: Path, assembly: str, fasta: Path | None, fai: Path | None, projected: Path | None, manual: Path | None, merged: Path | None, cases: list[dict[str, str]]) -> None:
    def uri(path: Path | None) -> str | None:
        return path.name if path else None

    tracks = []
    for track_id, name, path in (
        ("projected", "Projected GENCODE", projected),
        ("manual", "HAVANA manual", manual),
        ("merged", "Integrated annotation", merged),
    ):
        if path:
            tracks.append({
                "type": "FeatureTrack", "trackId": track_id, "name": name,
                "assemblyNames": [assembly],
                "adapter": {"type": "Gff3Adapter", "gffLocation": {"uri": uri(path)}},
            })
    displayed = []
    for case in cases:
        displayed.append({
            "refName": case["seqname"], "start": max(0, int(case["start"]) - 500),
            "end": int(case["end"]) + 500, "reversed": case["strand"] == "-",
            "assemblyName": assembly,
        })
    config = {
        "assemblies": [{
            "name": assembly, "sequence": {
                "type": "ReferenceSequenceTrack", "trackId": f"{assembly}-reference",
                "adapter": {"type": "IndexedFastaAdapter", "fastaLocation": {"uri": uri(fasta)}, "faiLocation": {"uri": uri(fai)}}
            }
        }] if fasta and fai else [],
        "tracks": tracks,
    }
    session = {
        "formatVersion": 1,
        "session": {
            "name": "Manual versus projected annotation review",
            "views": [{
                "type": "LinearGenomeView", "displayedRegions": displayed,
                "tracks": [{"type": "FeatureTrack", "configuration": x["trackId"]} for x in tracks]
            }]
        }
    }
    output.with_name("jbrowse_config.json").write_text(json.dumps(config, indent=2) + "\n")
    output.write_text(json.dumps(session, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--projected-gff", type=Path)
    parser.add_argument("--manual-gff", type=Path)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--decision-tsv", type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--fai", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    projected_genes: list[Model] = []
    projected_tx: dict[str, Model] = {}
    if args.projected_gff:
        _, pg, projected_tx = read_gff(args.projected_gff, retain_records=False)
        projected_genes = [Model(gid, gid, r.seqname, r.start, r.end, r.strand, r.feature, [r]) for gid, rs in pg.items() for r in rs]
    projected_gene_index = build_overlap_index(projected_genes)
    projected_tx_index = build_overlap_index(projected_tx.values())
    manual_genes: dict[str, list[Record]] = {}
    manual_tx: dict[str, Model] = {}
    if args.manual_gff:
        _, manual_genes, manual_tx = read_gff(args.manual_gff)
    decisions = load_decisions(args.decision_tsv)
    rows: list[dict[str, str]] = []
    for model in manual_tx.values():
        candidate_genes = overlapping_models(model, projected_gene_index)
        candidate_tx = overlapping_models(model, projected_tx_index)
        category, gene_ids, tx_ids = classify(model, candidate_genes, candidate_tx)
        review_required, review_reason, recommended_action = decision_annotation(category, gene_ids, tx_ids)
        rows.append({
            "case_id": f"{model.gene_id}:{model.ident}", "manual_gene_id": model.gene_id,
            "manual_transcript_id": model.ident, "seqname": model.seqname,
            "start": str(model.start), "end": str(model.end), "strand": model.strand,
            "manual_feature": model.feature, "classification": category,
            "projected_gene_ids": ",".join(gene_ids), "projected_transcript_ids": ",".join(tx_ids),
            "action": decisions.get(model.ident, {
                "exact_exon_structure": "add_manual_preserved",
                "readthrough_excluded": "add_manual_preserved",
                "novel_locus": "add_manual_preserved",
            }.get(category, "add_manual_preserved")),
            "best_projected_explanation": "",
            "review_required": review_required,
            "review_reason": review_reason,
            "recommended_action": recommended_action,
            "best_projected_transcript_id": "",
            "manual_exon_count": str(len(model.exons)),
            "projected_exon_count": "",
            "manual_transcript_start": str(model.start),
            "manual_transcript_end": str(model.end),
            "projected_transcript_start": "",
            "projected_transcript_end": "",
            "manual_five_prime_utr_bp": str(interval_span(model.five_utr)),
            "projected_five_prime_utr_bp": "",
            "manual_three_prime_utr_bp": str(interval_span(model.three_utr)),
            "projected_three_prime_utr_bp": "",
        })
        candidate_models = [candidate for candidate in candidate_tx if candidate.ident in tx_ids]
        if candidate_models:
            def match_rank(candidate: Model) -> tuple[int, int]:
                overlap = max(0, min(model.end, candidate.end) - max(model.start, candidate.start) + 1)
                if model.exons and candidate.exons and model.exons == candidate.exons:
                    rank = 3
                elif (model.exons and candidate.exons and len(model.exons) == len(candidate.exons)
                      and model.exons[1:-1] == candidate.exons[1:-1]):
                    rank = 2
                elif overlap / max(1, model.end - model.start + 1) >= 0.8:
                    rank = 1
                else:
                    rank = 0
                return rank, -abs(model.start - candidate.start) - abs(model.end - candidate.end)
            best = max(candidate_models, key=match_rank)
            rows[-1]["best_projected_explanation"] = model_explanation(model, best)
            rows[-1].update(model_comparison(model, best))
    review = args.output_dir / "annotation_review.tsv"
    with review.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["case_id"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    cases = args.output_dir / "review_cases.csv"
    with cases.open("w", newline="") as handle:
        fields = ["case_id", "label", "classification", "seqname", "start", "end", "strand", "manual_transcript_id", "projected_transcript_ids"]
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for row in rows:
            writer.writerow({"case_id": row["case_id"], "label": row["manual_transcript_id"], **{k: row[k] for k in fields[2:]}})
    representative = []
    for category in ("novel_locus", "same_locus_different_exon_chain", "overlapping_different_exon_chain", "same_internal_exon_chain_UTR_difference", "exact_exon_structure"):
        representative.extend([row for row in rows if row["classification"] == category][:3])
    representative_path = args.output_dir / "review_cases_representative.csv"
    with representative_path.open("w", newline="") as handle:
        fields = ["case_id", "label", "classification", "seqname", "start", "end", "strand", "manual_transcript_id", "projected_transcript_ids"]
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for row in representative:
            writer.writerow({"case_id": row["case_id"], "label": row["manual_transcript_id"], **{k: row[k] for k in fields[2:]}})
    handback = args.output_dir / "havana_decisions.tsv"
    handback_fields = [
        "case_id", "manual_gene_id", "manual_transcript_id", "seqname", "start", "end", "strand",
        "classification", "projected_gene_ids", "projected_transcript_ids", "best_projected_explanation",
        "review_reason", "recommended_action", "havana_decision", "havana_comment",
    ]
    with handback.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=handback_fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            if row["review_required"] != "true":
                continue
            writer.writerow({**{field: row.get(field, "") for field in handback_fields[:-2]}, "havana_decision": "", "havana_comment": ""})
    merged = args.output_dir / "integrated.gff3"
    drop_projected = {
        tx_id
        for row in rows if row["action"] == "keep_manual"
        for tx_id in row["projected_transcript_ids"].split(",") if tx_id
    }
    curated_projected: dict[str, list[str]] = defaultdict(list)
    manual_metadata: dict[str, dict[str, str]] = {}
    gene_metadata: dict[str, dict[str, str]] = {}
    effective_decisions = dict(decisions)
    for row in rows:
        effective_decisions[row["manual_transcript_id"]] = row["action"]
        manual_metadata[row["manual_transcript_id"]] = {
            "manual_review_required": row["review_required"],
            "manual_review_reason": row["review_reason"],
            "manual_recommended_action": row["recommended_action"],
            "manual_policy_classification": row["classification"],
        }
        gene_metadata.setdefault(row["manual_gene_id"], {
            "manual_review_required": row["review_required"],
            "manual_review_reason": row["review_reason"],
            "manual_recommended_action": row["recommended_action"],
            "manual_policy_classification": row["classification"],
        })
        if row["action"] == "mark_manual_curated":
            for tx_id in row["projected_transcript_ids"].split(","):
                if tx_id:
                    curated_projected[tx_id].append(row["manual_transcript_id"])
            for gene_id in row["projected_gene_ids"].split(","):
                if gene_id:
                    curated_projected[gene_id].append(row["manual_transcript_id"])
    write_merged(args.projected_gff, args.manual_gff, manual_genes, manual_tx, effective_decisions, drop_projected, curated_projected, manual_metadata, gene_metadata, merged)
    write_session(args.output_dir / "jbrowse_session.json", args.assembly, args.fasta, args.fai, args.projected_gff, args.manual_gff, merged, representative)
    (args.output_dir / "README.txt").write_text("Use review_cases.csv to select cases. jbrowse_config.json registers the assembly and tracks; jbrowse_session.json opens the first review cases.\n")


if __name__ == "__main__":
    main()
