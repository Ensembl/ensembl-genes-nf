#!/usr/bin/env python3
"""Build an auditable translon_db source manifest from TransCODE pilot outputs.

The manifest is deliberately provenance-first.  It starts from the canonical
pilot sample space, maps raw tool leaves onto canonical keys, records native
caller classes without interpreting them as CDS/non-CDS truth, and preserves
non-ingested leaves as explicit excluded rows.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


TOOLS = ("iRibo", "ORFQuant", "RibORF2", "RiboTIE", "PRICE")
FASTQ_ROUTE_TOOLS = ("iRibo", "ORFQuant", "RibORF2", "RiboTIE", "PRICE")
NATIVE_CLASSES = ("annotated", "novel")
MATCHED = "matched"
RIBORF2_FASTQ_DIRNAME = "fastq_RibORF2.0_ORFidentification"


@dataclass(frozen=True)
class SampleMeta:
    sample_no: int
    sample_name: str
    dataset_id: str
    sample_id: str
    cell_type: str
    replicate: str
    is_pooled: bool
    run_type: str
    fastq_route_expected: int


EXPECTED_SAMPLE_METADATA: tuple[SampleMeta, ...] = (
    SampleMeta(1, "Pancreas 1", "GSE144682", "SRR11005875_to_79", "pancreas", "1", False, "individual_sample", 1),
    SampleMeta(2, "Pancreas 2", "GSE144682", "SRR11005880_to_84", "pancreas", "2", False, "individual_sample", 1),
    SampleMeta(3, "Pancreas 3", "GSE144682", "SRR11005885_to_89", "pancreas", "3", False, "individual_sample", 1),
    SampleMeta(4, "Pancreas 4", "GSE144682", "SRR11005890_to_94", "pancreas", "4", False, "individual_sample", 1),
    SampleMeta(5, "Pancreas 5", "GSE144682", "SRR11005895_to_99", "pancreas", "5", False, "individual_sample", 1),
    SampleMeta(6, "Pancreas 6", "GSE144682", "SRR11005900_to_04", "pancreas", "6", False, "individual_sample", 1),
    SampleMeta(7, "Fibroblast 1", "GSE182371", "SRR15513179", "fibroblast", "1", False, "individual_sample", 0),
    SampleMeta(8, "Fibroblast 2", "GSE182371", "SRR15513180", "fibroblast", "2", False, "individual_sample", 0),
    SampleMeta(9, "Fibroblast 3", "GSE182371", "SRR15513181", "fibroblast", "3", False, "individual_sample", 0),
    SampleMeta(10, "Fibroblast 4", "GSE182371", "SRR15513182", "fibroblast", "4", False, "individual_sample", 0),
    SampleMeta(11, "Fibroblast 5", "GSE182371", "Fib_24_45m", "fibroblast", "5", False, "individual_sample", 0),
    SampleMeta(12, "Fibroblast 6", "GSE182371", "Fib_24_bsl", "fibroblast", "6", False, "individual_sample", 0),
    SampleMeta(13, "Fibroblast 7", "GSE182371", "Fib_27_45m", "fibroblast", "7", False, "individual_sample", 0),
    SampleMeta(14, "Fibroblast 8", "GSE182371", "Fib_27_bsl", "fibroblast", "8", False, "individual_sample", 0),
    SampleMeta(15, "Fibroblast 9", "GSE182371", "Fib_41_45m", "fibroblast", "9", False, "individual_sample", 0),
    SampleMeta(16, "Fibroblast 10", "GSE182371", "Fib_41_bsl", "fibroblast", "10", False, "individual_sample", 0),
    SampleMeta(17, "Endothelial cell 1", "GSE182371", "SRR15513197", "endothelial", "1", False, "individual_sample", 0),
    SampleMeta(18, "Endothelial cell 2", "GSE182371", "SRR15513198_GENELAB1026", "endothelial", "2", False, "individual_sample", 0),
    SampleMeta(19, "Endothelial cell 3", "GSE182371", "SRR15513199", "endothelial", "3", False, "individual_sample", 0),
    SampleMeta(20, "Endothelial cell 4", "GSE182371", "SRR15513200", "endothelial", "4", False, "individual_sample", 0),
    SampleMeta(21, "Endothelial cell 5", "GSE182371", "SRR15513201", "endothelial", "5", False, "individual_sample", 0),
    SampleMeta(22, "Endothelial cell 6", "GSE182371", "SRR15513202_GENELAB1143", "endothelial", "6", False, "individual_sample", 0),
    SampleMeta(23, "Pooled Pancreas cells", "GSE144682", "Ribo_Pancreas_pooled", "pancreas", "pooled", True, "pooled_or_aggregate", 0),
    SampleMeta(24, "Pooled Fibroblast cells", "GSE182371", "Ribo_Fib_pooled", "fibroblast", "pooled", True, "pooled_or_aggregate", 0),
    SampleMeta(25, "Pooled Endothelial cells", "GSE182371", "Ribo_EC_pooled", "endothelial", "pooled", True, "pooled_or_aggregate", 0),
)
EXPECTED_SAMPLE_BY_ID = {sample.sample_id: sample for sample in EXPECTED_SAMPLE_METADATA}
EXPECTED_SAMPLES = frozenset(EXPECTED_SAMPLE_BY_ID)


PRICE_SAMPLE_MAP: dict[str, str] = {
    "Fibo": "Ribo_Fib_pooled",
    "Fibo_0": "SRR15513179",
    "Fibo_1": "SRR15513180",
    "Fibo_2": "SRR15513181",
    "Fibo_3": "SRR15513182",
    "Fibo_4": "Fib_24_45m",
    "Fibo_5": "Fib_24_bsl",
    "Fibo_6": "Fib_27_45m",
    "Fibo_7": "Fib_27_bsl",
    "Fibo_8": "Fib_41_45m",
    "Fibo_9": "Fib_41_bsl",
    "Endo": "Ribo_EC_pooled",
    "Endo_0": "SRR15513197",
    "Endo_1": "SRR15513198_GENELAB1026",
    "Endo_2": "SRR15513199",
    "Endo_3": "SRR15513200",
    "Endo_4": "SRR15513201",
    "Endo_5": "SRR15513202_GENELAB1143",
    "Pancreas": "Ribo_Pancreas_pooled",
    "Pancreas_0": "SRR11005875_to_79",
    "Pancreas_1": "SRR11005880_to_84",
    "Pancreas_2": "SRR11005885_to_89",
    "Pancreas_3": "SRR11005890_to_94",
    "Pancreas_4": "SRR11005895_to_99",
    "Pancreas_5": "SRR11005900_to_04",
    "pancreas_mymapping": "Ribo_Pancreas_pooled",
    "pancreas_mymapping_0": "SRR11005875_to_79",
    "pancreas_mymapping_1": "SRR11005880_to_84",
    "pancreas_mymapping_2": "SRR11005885_to_89",
    "pancreas_mymapping_3": "SRR11005890_to_94",
    "pancreas_mymapping_4": "SRR11005895_to_99",
    "pancreas_mymapping_5": "SRR11005900_to_04",
}
IRIBO_POOLED_MAP = {
    "Ribo_EC_p": "Ribo_EC_pooled",
    "Ribo_Fib_p": "Ribo_Fib_pooled",
    "Ribo_pancreas_p": "Ribo_Pancreas_pooled",
}
MISC_FIXES = {
    "Ribo_pancreas_pooled": "Ribo_Pancreas_pooled",
    "Ribo_ECs_pooled": "Ribo_EC_pooled",
    "Ribo_fib_pooled": "Ribo_Fib_pooled",
}
RIBORF2_FASTQ_SAMPLE_MAP = {
    "Pancreas1": "SRR11005875_to_79",
    "Pancreas2": "SRR11005880_to_84",
    "Pancreas3": "SRR11005885_to_89",
    "Pancreas4": "SRR11005890_to_94",
    "Pancreas5": "SRR11005895_to_99",
    "Pancreas6": "SRR11005900_to_04",
    "Pooledalignmentfile": "Ribo_Pancreas_pooled",
}


FIELDNAMES = [
    "sample_id",
    "sample_no",
    "sample_name",
    "dataset_id",
    "cell_type",
    "replicate",
    "is_pooled",
    "run_type",
    "tool",
    "route",
    "input_route",
    "native_class",
    "source_feature_class",
    "path",
    "source_path",
    "raw_record_count",
    "ingest_status",
    "reason",
    "reconciliation",
    "parser",
    "fastq_route_expected",
]


def normalise(raw: str) -> str:
    sample = raw.strip()
    sample = sample.replace("GENELAB-000", "GENELAB").replace("GENELAB_000", "GENELAB")
    sample = PRICE_SAMPLE_MAP.get(sample, sample)
    sample = IRIBO_POOLED_MAP.get(sample, sample)
    sample = MISC_FIXES.get(sample, sample)
    sample = re.sub(r"_1$", "", sample)
    return sample


def route_to_input_route(route: str) -> str:
    return "fastq_to_orf" if route == "fastq" else "bam_to_orf"


def class_to_source_feature_class(native_class: str) -> str:
    # Tool-native annotated/known/novel classes are provenance labels, not the
    # central GENCODE-derived CDS/non-CDS assignment.
    return "unknown"


def parser_for(tool: str, path: Path) -> str:
    if tool == "iRibo":
        return "iribo_gff"
    if tool == "ORFQuant":
        return "orfquant_bed_exon"
    if tool == "RibORF2" or path.suffix == ".bed12":
        return "bed12"
    if tool == "RiboTIE":
        return "translonscorer_csv" if path.suffix == ".csv" else "ribotie_gtf"
    if tool == "PRICE":
        return "bed12"
    return "unknown"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def count_records(path: Path) -> int:
    """Count data records, excluding comments, browser tracks, and CSV headers."""
    try:
        with open_text(path) as handle:
            count = 0
            first_data = True
            for line in handle:
                if not line.strip() or line.startswith("#") or line.startswith("track"):
                    continue
                if path.suffix == ".csv" and first_data:
                    first_data = False
                    continue
                first_data = False
                count += 1
    except OSError:
        return 0
    return count


def with_metadata(row: dict[str, object]) -> dict[str, object]:
    sample_id = str(row.get("sample_id", ""))
    meta = EXPECTED_SAMPLE_BY_ID.get(sample_id)
    base = {field: "" for field in FIELDNAMES}
    base.update(row)
    base["source_path"] = base.get("source_path") or base.get("path", "")
    base["input_route"] = base.get("input_route") or route_to_input_route(str(base.get("route", "bam")))
    base["source_feature_class"] = base.get("source_feature_class") or class_to_source_feature_class(str(base.get("native_class", "")))
    if meta:
        base.update(
            {
                "sample_no": meta.sample_no,
                "sample_name": meta.sample_name,
                "dataset_id": meta.dataset_id,
                "cell_type": meta.cell_type,
                "replicate": meta.replicate,
                "is_pooled": int(meta.is_pooled),
                "run_type": meta.run_type,
                "fastq_route_expected": meta.fastq_route_expected,
            }
        )
    return base


def matched_row(tool: str, sample_id: str, route: str, native_class: str, path: Path) -> dict[str, object]:
    return with_metadata(
        {
            "sample_id": sample_id,
            "tool": tool,
            "route": route,
            "native_class": native_class,
            "path": str(path),
            "source_path": str(path),
            "raw_record_count": count_records(path),
            "ingest_status": MATCHED,
            "reason": "",
            "parser": parser_for(tool, path),
        }
    )


def status_row(
    tool: str,
    sample_id: str,
    route: str,
    native_class: str,
    ingest_status: str,
    reason: str,
) -> dict[str, object]:
    return with_metadata(
        {
            "sample_id": sample_id,
            "tool": tool,
            "route": route,
            "native_class": native_class,
            "ingest_status": ingest_status,
            "reason": reason,
            "raw_record_count": "",
        }
    )


def excluded_row(path: Path, reason: str) -> dict[str, object]:
    return with_metadata(
        {
            "path": str(path),
            "source_path": str(path),
            "ingest_status": "excluded",
            "reason": reason,
        }
    )


def collect_iribo(results_dir: Path) -> tuple[list[dict[str, object]], set[Path]]:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()
    for path in sorted((results_dir / "iRibo_results").glob("*.bed")):
        match = re.match(r"^(annotated|unannotated)_orfs_(.+?)(?:_S\d+_R1_001)?\.bed$", path.name)
        if not match:
            continue
        prefix, raw_sample = match.groups()
        route = "fastq" if raw_sample.endswith("_fastq") else "bam"
        raw_sample = re.sub(r"_fastq$", "", raw_sample)
        native_class = "annotated" if prefix == "annotated" else "novel"
        rows.append(matched_row("iRibo", normalise(raw_sample), route, native_class, path))
        seen.add(path)
    return rows, seen


def collect_orfquant(results_dir: Path) -> tuple[list[dict[str, object]], set[Path]]:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()
    for path in sorted((results_dir / "ORFQuant_results" / "bed_files").glob("*.bed")):
        match = re.match(
            r"^(annotated|novel)_orf_exon_genomic_"
            r"(.+?)(?:_S\d+_R1_001(?:_trimmed)?)?(?:\.Aligned\.sortedByCoord\.out)?\.bed$",
            path.name,
        )
        if not match:
            continue
        cls, raw_sample = match.groups()
        native_class = "annotated" if cls == "annotated" else "novel"
        sample_id = normalise(raw_sample)
        if ".Aligned.sortedByCoord.out.bed" in path.name or "_trimmed.Aligned.sortedByCoord.out.bed" in path.name:
            route = "bam"
        elif sample_id in {sample.sample_id for sample in EXPECTED_SAMPLE_METADATA[:6]}:
            route = "fastq"
        else:
            # Pooled ORFQuant outputs do not carry STAR's BAM stem, but the
            # staged input design only has pooled BAMs, not pooled FASTQs.
            route = "bam"
        rows.append(matched_row("ORFQuant", sample_id, route, native_class, path))
        seen.add(path)
    for path in sorted((results_dir / "ORFQuant_results" / "bed_files").glob("*.bed")):
        if path not in seen and re.match(r"^(annotated|novel)_orf_exon_genomic_", path.name):
            rows.append(excluded_row(path, "orfquant_noncanonical_bed"))
            seen.add(path)
    return rows, seen


def collect_riborf2(results_dir: Path, converted_dir: Path) -> tuple[list[dict[str, object]], set[Path]]:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()
    if converted_dir.exists():
        for path in sorted(converted_dir.glob("*.bed12")):
            stem = re.sub(r"_S\d+_R1_001(?:_trimmed)?$", "", path.stem)
            rows.append(matched_row("RibORF2", normalise(stem), "bam", "unknown", path))
            seen.add(path)

    matched_fastq_keys: set[tuple[str, str]] = set()
    fastq_result_dir = results_dir / RIBORF2_FASTQ_DIRNAME
    if fastq_result_dir.exists():
        for path in sorted(fastq_result_dir.rglob("*.bed")):
            sample_dir = path.parent.name
            sample_id = RIBORF2_FASTQ_SAMPLE_MAP.get(sample_dir)
            if not sample_id:
                continue
            lower_name = path.name.lower()
            if "annotatedorf" in lower_name:
                native_class = "annotated"
            elif "novelsmorf" in lower_name or "novelsmorfs" in lower_name:
                native_class = "novel"
            else:
                continue
            rows.append(matched_row("RibORF2", sample_id, "fastq", native_class, path))
            seen.add(path)
            matched_fastq_keys.add((sample_id, native_class))

    legacy_fastq_dir = results_dir / "RibORF_results" / "RibORF_Output_fastqtoORFcalling"
    if not matched_fastq_keys and (not legacy_fastq_dir.exists() or not any(legacy_fastq_dir.rglob("*"))):
        for sample in EXPECTED_SAMPLE_METADATA[:6]:
            for native_class in NATIVE_CLASSES:
                rows.append(status_row("RibORF2", sample.sample_id, "fastq", native_class, "empty_at_source", "RibORF2 FASTQ output directory is empty"))
    return rows, seen


def ribotie_sample_and_route(stem: str) -> tuple[str, str]:
    stem = re.sub(r"^db_", "", stem)
    stem = re.sub(r"\.(annotated|novel)\.out$", "", stem)
    stem = re.sub(r"\.Aligned\.toTranscriptome\.out$", "", stem)
    stem = re.sub(r"_Transcriptome$", "", stem)
    pancreas = re.match(r"^Pancreas_([1-6])$", stem)
    if pancreas:
        idx = int(pancreas.group(1)) - 1
        return EXPECTED_SAMPLE_METADATA[idx].sample_id, "fastq"
    return normalise(stem), "bam"


def collect_ribotie(results_dir: Path) -> tuple[list[dict[str, object]], set[Path]]:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()
    root = results_dir / "RiboTIE_results" / "deliverables"
    candidates: dict[tuple[str, str, str], list[tuple[int, Path]]] = {}

    def add_candidate(path: Path, native_class: str, priority: int) -> None:
        if "unfiltered" in path.as_posix():
            return
        sample_id, route = ribotie_sample_and_route(path.stem)
        if route == "fastq":
            rows.append(excluded_row(path, "ribotie_fastq_split_deliverable_superseded_by_raw_fastq_gtf"))
            seen.add(path)
            return
        candidates.setdefault((sample_id, route, native_class), []).append((priority, path))

    for subdir, native_class in (("annotated", "annotated"), ("novel", "novel")):
        for path in sorted((root / subdir).glob("*.gtf")):
            # Prefer the delivered .Aligned.toTranscriptome.out GTFs.  Despite
            # the name, these carry parseable RiboTIE CDS blocks in the current
            # pilot delivery; the shorter db_SAMPLE.annotated.out.gtf files can
            # be empty/non-ORF after conversion for some samples.
            priority = 0 if ".Aligned.toTranscriptome.out." in path.name else 1
            add_candidate(path, native_class, priority)
        for path in sorted((root / subdir).glob("*.csv")):
            add_candidate(path, native_class, 2)

    raw_fastq = results_dir / "RiboTIE_results" / "raw" / "fastq"
    for path in sorted(raw_fastq.glob("*.gtf")):
        if "unfiltered" in path.name:
            continue
        sample_id, route = ribotie_sample_and_route(path.stem)
        rows.append(matched_row("RiboTIE", sample_id, route, "unknown", path))
        seen.add(path)
    for path in sorted(raw_fastq.glob("*.csv")):
        if "unfiltered" in path.name:
            continue
        rows.append(excluded_row(path, "ribotie_fastq_csv_superseded_by_gtf"))
        seen.add(path)

    for key, ranked_paths in sorted(candidates.items()):
        ranked_paths = sorted(ranked_paths, key=lambda item: (item[0], item[1].name))
        best_priority, best_path = ranked_paths[0]
        rows.append(matched_row("RiboTIE", key[0], key[1], key[2], best_path))
        seen.add(best_path)
        for priority, path in ranked_paths[1:]:
            if priority == best_priority:
                reason = "ribotie_duplicate_equivalent_deliverable"
            elif path.suffix == ".csv":
                reason = "ribotie_csv_superseded_by_split_gtf"
            else:
                reason = "ribotie_nonpreferred_gtf_superseded_by_aligned_gtf"
            rows.append(excluded_row(path, reason))
            seen.add(path)

    for path in sorted(root.rglob("*")):
        if path.is_file() and path not in seen:
            rows.append(excluded_row(path, "ribotie_noncanonical_deliverable_or_unfiltered"))
            seen.add(path)
    return rows, seen


def collect_price(results_dir: Path) -> tuple[list[dict[str, object]], set[Path]]:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()
    for path in sorted((results_dir / "PRICE_results" / "price").glob("*.bed")):
        match = re.match(r"^(.+)\.(known|novel)\.bed$", path.name)
        if not match:
            continue
        raw_sample, cls = match.groups()
        route = "fastq" if raw_sample.startswith("pancreas_mymapping") else "bam"
        native_class = "annotated" if cls == "known" else "novel"
        rows.append(matched_row("PRICE", normalise(raw_sample), route, native_class, path))
        seen.add(path)
    return rows, seen


def iter_all_leaves(results_dir: Path, converted_dir: Path) -> Iterable[Path]:
    if results_dir.exists():
        yield from (path for path in results_dir.rglob("*") if path.is_file())
    if converted_dir.exists() and results_dir not in converted_dir.parents and converted_dir != results_dir:
        yield from (path for path in converted_dir.rglob("*") if path.is_file())


def design_matrix_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample in EXPECTED_SAMPLE_METADATA:
        for tool in TOOLS:
            native_classes = ("unknown",) if tool == "RibORF2" else NATIVE_CLASSES
            for native_class in native_classes:
                rows.append(status_row(tool, sample.sample_id, "bam", native_class, "expected", ""))
        if sample.fastq_route_expected:
            for tool in FASTQ_ROUTE_TOOLS:
                native_classes = ("unknown",) if tool == "RiboTIE" else NATIVE_CLASSES
                for native_class in native_classes:
                    rows.append(status_row(tool, sample.sample_id, "fastq", native_class, "expected", ""))
    return rows


def optional_key(key: tuple[str, str, str, str]) -> bool:
    _tool, sample_id, route, _native_class = key
    return sample_id == "Ribo_Pancreas_pooled" and route == "fastq"


def is_allowed_missing_cell(key: tuple[str, str, str, str], ingest_status: str) -> bool:
    tool, sample_id, route, _native_class = key
    if route == "fastq" and sample_id not in {sample.sample_id for sample in EXPECTED_SAMPLE_METADATA[:6]}:
        return ingest_status == "empty_at_source"
    if route == "fastq" and tool == "RibORF2":
        return ingest_status == "empty_at_source"
    return False


def expected_key_set() -> set[tuple[str, str, str, str]]:
    return {
        (str(row["tool"]), str(row["sample_id"]), str(row["route"]), str(row["native_class"]))
        for row in design_matrix_rows()
    }


def add_expected_missing(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    present = {
        (str(row["tool"]), str(row["sample_id"]), str(row["route"]), str(row["native_class"]))
        for row in rows
        if row.get("ingest_status") == MATCHED
    }
    existing_status = {
        (str(row["tool"]), str(row["sample_id"]), str(row["route"]), str(row["native_class"]))
        for row in rows
        if row.get("ingest_status") in {"empty_at_source", "blocked"}
    }
    for expected in design_matrix_rows():
        key = (
            str(expected["tool"]),
            str(expected["sample_id"]),
            str(expected["route"]),
            str(expected["native_class"]),
        )
        if key in present or key in existing_status:
            continue
        status = "blocked"
        reason = "expected_by_design_but_no_source_file_matched"
        if expected["route"] == "fastq" and expected["sample_id"] not in {
            *(sample.sample_id for sample in EXPECTED_SAMPLE_METADATA[:6]),
            "Ribo_Pancreas_pooled",
        }:
            status = "empty_at_source"
            reason = "FASTQ route was not requested for this sample"
        rows.append(
            status_row(
                str(expected["tool"]),
                str(expected["sample_id"]),
                str(expected["route"]),
                str(expected["native_class"]),
                status,
                reason,
            )
        )
    return rows


def reconcile(rows: list[dict[str, object]], old_manifest: Path | None, parser_manifest: Path | None) -> None:
    old_keys: set[tuple[str, str, str, str]] = set()
    for manifest in (old_manifest, parser_manifest):
        if not manifest or not manifest.exists():
            continue
        opener = gzip.open if manifest.suffix == ".gz" else open
        with opener(manifest, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                tool = row.get("tool") or row.get("tool_hint") or row.get("detected_tool") or ""
                sample_id = row.get("sample_id", "")
                route = row.get("route") or ("fastq" if row.get("input_route") == "fastq_to_orf" else "bam")
                native_class = row.get("native_class") or ""
                if tool and sample_id:
                    old_keys.add((tool, sample_id, route, native_class))
    if not old_keys:
        for row in rows:
            row["reconciliation"] = ""
        return
    for row in rows:
        if row.get("ingest_status") != MATCHED:
            row["reconciliation"] = ""
            continue
        key = (str(row["tool"]), str(row["sample_id"]), str(row["route"]), str(row["native_class"]))
        loose_key = (str(row["tool"]), str(row["sample_id"]), str(row["route"]), "")
        row["reconciliation"] = "agree" if key in old_keys or loose_key in old_keys else "manifest_missing"


def assert_manifest(rows: list[dict[str, object]], matched_paths: set[Path], all_leaves: set[Path]) -> None:
    errors: list[str] = []
    expected_keys = expected_key_set()
    keyed_rows = [
        (row["tool"], row["sample_id"], row["route"], row["native_class"])
        for row in rows
        if row.get("tool") and row.get("sample_id") and row.get("route") and row.get("native_class")
    ]
    matched_keys = [
        (row["tool"], row["sample_id"], row["route"], row["native_class"])
        for row in rows
        if row.get("ingest_status") == MATCHED
    ]
    unexpected_keys = sorted(key for key in set(keyed_rows) - expected_keys if not optional_key(key))
    if unexpected_keys:
        errors.append(f"rows outside expected design matrix: {unexpected_keys[:10]}")
    missing_keys = sorted(expected_keys - set(keyed_rows))
    if missing_keys:
        errors.append(f"expected design cells absent from manifest: {missing_keys[:10]}")
    duplicates = [key for key, count in Counter(matched_keys).items() if count > 1]
    if duplicates:
        errors.append(f"duplicate matched keys: {duplicates[:10]}")
    status_by_key = {
        (row["tool"], row["sample_id"], row["route"], row["native_class"]): row
        for row in rows
        if row.get("tool") and row.get("sample_id") and row.get("route") and row.get("native_class")
    }
    unexpected_missing = [
        key for key, row in status_by_key.items()
        if row.get("ingest_status") in {"blocked", "empty_at_source"}
        and not is_allowed_missing_cell(key, str(row.get("ingest_status", "")))
    ]
    if unexpected_missing:
        errors.append(f"unexpected missing source cells: {unexpected_missing[:10]}")
    bad_status_rows = [
        row for row in rows
        if row.get("ingest_status") not in {MATCHED, "blocked", "empty_at_source", "excluded"}
    ]
    if bad_status_rows:
        errors.append(f"rows with invalid ingest_status: {len(bad_status_rows)}")
    missing_reasons = [
        row for row in rows
        if row.get("ingest_status") != MATCHED and not row.get("reason")
    ]
    if missing_reasons:
        errors.append(f"non-matched rows missing reason: {len(missing_reasons)}")
    matched_missing_path = [
        row for row in rows
        if row.get("ingest_status") == MATCHED and (not row.get("source_path") or not Path(str(row["source_path"])).exists())
    ]
    if matched_missing_path:
        errors.append(f"matched rows with missing source_path on disk: {len(matched_missing_path)}")
    matched_without_records = [
        row for row in rows
        if row.get("ingest_status") == MATCHED and str(row.get("raw_record_count", "")) in {"", "0"}
    ]
    if matched_without_records:
        errors.append(f"matched rows with zero raw_record_count: {len(matched_without_records)}")
    matched_not_unknown_class = [
        row for row in rows
        if row.get("ingest_status") == MATCHED and row.get("source_feature_class") != "unknown"
    ]
    if matched_not_unknown_class:
        errors.append(f"matched rows should keep source_feature_class=unknown before DB classification: {len(matched_not_unknown_class)}")
    unmapped_matched = [row for row in rows if row.get("ingest_status") == MATCHED and row.get("sample_id") not in EXPECTED_SAMPLES]
    if unmapped_matched:
        errors.append(f"matched rows with noncanonical sample_id: {len(unmapped_matched)}")
    illegal_fastq = [
        row for row in rows
        if row.get("ingest_status") == MATCHED
        and row.get("route") == "fastq"
        and row.get("sample_id") not in {*(sample.sample_id for sample in EXPECTED_SAMPLE_METADATA[:6]), "Ribo_Pancreas_pooled"}
    ]
    if illegal_fastq:
        errors.append(f"illegal FASTQ matched rows: {len(illegal_fastq)}")
    price_zero_index = {
        ("PRICE", "SRR11005875_to_79", "bam", "annotated"),
        ("PRICE", "SRR11005875_to_79", "fastq", "annotated"),
    }
    if not price_zero_index.issubset(set(matched_keys)):
        errors.append("PRICE zero-index sanity failed for Pancreas_0/pancreas_mymapping_0")
    accounted = {
        Path(str(row["source_path"]))
        for row in rows
        if row.get("source_path") and row.get("ingest_status") in {MATCHED, "excluded"}
    }
    missing_accounting = all_leaves - accounted
    if missing_accounting:
        errors.append(f"unaccounted leaves under known output roots: {len(missing_accounting)}")
    if errors:
        raise SystemExit("Manifest assertions failed:\n  - " + "\n  - ".join(errors))


def print_summary(rows: list[dict[str, object]]) -> None:
    counts = Counter(str(row["ingest_status"]) for row in rows)
    matched_by_tool = Counter(str(row["tool"]) for row in rows if row.get("ingest_status") == MATCHED)
    print("\nManifest rows by status:")
    for status, count in sorted(counts.items()):
        print(f"  {status}: {count}")
    print("Matched rows by tool:")
    for tool in TOOLS:
        print(f"  {tool}: {matched_by_tool[tool]}")
    blocked = [row for row in rows if row.get("ingest_status") in {"blocked", "empty_at_source"}]
    print(f"Expected cells without matched source: {len(blocked)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--riborf2-converted", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--design-matrix-out", type=Path)
    parser.add_argument("--old-manifest", type=Path)
    parser.add_argument("--parser-manifest", type=Path)
    parser.add_argument("--no-assert", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    collectors = (
        collect_iribo(args.results_dir),
        collect_orfquant(args.results_dir),
        collect_riborf2(args.results_dir, args.riborf2_converted),
        collect_ribotie(args.results_dir),
        collect_price(args.results_dir),
    )
    rows: list[dict[str, object]] = []
    matched_or_excluded_paths: set[Path] = set()
    for tool_rows, seen in collectors:
        rows.extend(tool_rows)
        matched_or_excluded_paths.update(seen)

    all_leaves = set(iter_all_leaves(args.results_dir, args.riborf2_converted))
    for path in sorted(all_leaves - matched_or_excluded_paths):
        rows.append(excluded_row(path, "outside_canonical_ingest_rules"))

    rows = add_expected_missing(rows)
    reconcile(rows, args.old_manifest, args.parser_manifest)
    if not args.no_assert:
        assert_manifest(rows, matched_or_excluded_paths, all_leaves)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    if args.design_matrix_out:
        args.design_matrix_out.parent.mkdir(parents=True, exist_ok=True)
        with args.design_matrix_out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDNAMES, lineterminator="\n")
            writer.writeheader()
            writer.writerows(design_matrix_rows())

    print(f"Wrote {len(rows)} rows -> {args.out}")
    if args.design_matrix_out:
        print(f"Wrote design matrix -> {args.design_matrix_out}")
    print_summary(rows)


if __name__ == "__main__":
    main()
