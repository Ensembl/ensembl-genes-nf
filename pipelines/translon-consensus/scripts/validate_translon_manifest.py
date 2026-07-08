#!/usr/bin/env python3
"""Stage 1 datachecks: validate a translon manifest against design and disk."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_translon_manifest import (  # noqa: E402
    EXPECTED_SAMPLE_METADATA,
    FASTQ_ROUTE_TOOLS,
    MATCHED,
    TOOLS,
    expected_key_set,
    is_allowed_missing_cell,
    iter_all_leaves,
    optional_key,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row.get("tool", ""), row.get("sample_id", ""), row.get("route", ""), row.get("native_class", ""))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--riborf2-converted", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--reconciliation-dispositions", type=Path)
    parser.add_argument("--ignore-prefix", type=Path, action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_tsv(args.manifest)
    errors: list[str] = []
    ignore_prefixes = [path.resolve() for path in args.ignore_prefix]

    def ignored(path: Path) -> bool:
        resolved = path.resolve()
        return any(resolved == prefix or prefix in resolved.parents for prefix in ignore_prefixes)

    manifest_paths = [
        Path(row["source_path"])
        for row in rows
        if row.get("source_path") and row.get("ingest_status") in {MATCHED, "excluded"}
        and not ignored(Path(row["source_path"]))
    ]
    path_counts = Counter(manifest_paths)
    duplicated_paths = [str(path) for path, count in path_counts.items() if count > 1]
    if duplicated_paths:
        errors.append(f"source paths appearing more than once as matched/excluded: {duplicated_paths[:20]}")
    missing_on_disk = [str(path) for path in manifest_paths if not path.exists()]
    if missing_on_disk:
        errors.append(f"manifest source_path values missing on disk: {missing_on_disk[:20]}")

    disk_leaves = {path for path in iter_all_leaves(args.results_dir, args.riborf2_converted) if not ignored(path)}
    manifest_leaf_set = set(manifest_paths)
    missing_from_manifest = sorted(str(path) for path in disk_leaves - manifest_leaf_set)
    phantom_paths = sorted(str(path) for path in manifest_leaf_set - disk_leaves)
    if missing_from_manifest:
        errors.append(f"disk leaves absent from manifest closure: {missing_from_manifest[:20]}")
    if phantom_paths:
        errors.append(f"manifest paths outside filesystem closure: {phantom_paths[:20]}")

    expected = expected_key_set()
    rows_by_key: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        k = key(row)
        if all(k):
            rows_by_key[k].append(row)
    unexpected = sorted(k for k in rows_by_key if k not in expected and not optional_key(k))
    if unexpected:
        errors.append(f"manifest rows outside design matrix: {unexpected[:20]}")
    absent_expected = sorted(k for k in expected if k not in rows_by_key)
    if absent_expected:
        errors.append(f"design cells absent from manifest: {absent_expected[:20]}")

    coverage_rows: list[dict[str, object]] = []
    all_grid_keys = sorted(expected | {k for k in rows_by_key if optional_key(k)})
    for k in all_grid_keys:
        cell_rows = rows_by_key.get(k, [])
        matched = [row for row in cell_rows if row.get("ingest_status") == MATCHED]
        missing = [row for row in cell_rows if row.get("ingest_status") in {"blocked", "empty_at_source"}]
        status = "matched" if matched else (missing[0].get("ingest_status", "absent") if missing else "absent")
        reason = "" if matched else (missing[0].get("reason", "") if missing else "")
        coverage_rows.append(
            {
                "tool": k[0],
                "sample_id": k[1],
                "route": k[2],
                "native_class": k[3],
                "status": status,
                "matched_rows": len(matched),
                "reason": reason,
            }
        )
        if status != "matched" and not is_allowed_missing_cell(k, status):
            errors.append(f"unexpected non-matched design cell: {k} status={status} reason={reason}")
        if status != "matched" and not reason:
            errors.append(f"non-matched design cell lacks reason: {k}")
    write_tsv(args.out_dir / "stage1_coverage_grid.tsv", coverage_rows, list(coverage_rows[0]))

    fastq_rows = [row for row in rows if row.get("ingest_status") == MATCHED and row.get("route") == "fastq"]
    allowed_fastq_samples = {sample.sample_id for sample in EXPECTED_SAMPLE_METADATA[:6]} | {"Ribo_Pancreas_pooled"}
    illegal_fastq = [row for row in fastq_rows if row.get("sample_id") not in allowed_fastq_samples]
    if illegal_fastq:
        errors.append(f"illegal FASTQ matched rows: {[row.get('source_path') for row in illegal_fastq[:20]]}")
    illegal_fastq_tools = [row for row in fastq_rows if row.get("tool") not in FASTQ_ROUTE_TOOLS]
    if illegal_fastq_tools:
        errors.append(f"FASTQ matched rows for tools without FASTQ route: {[row.get('source_path') for row in illegal_fastq_tools[:20]]}")

    for route in ("bam", "fastq"):
        price_key = ("PRICE", "SRR11005875_to_79", route, "annotated")
        paths = [row.get("source_path", "") for row in rows_by_key.get(price_key, []) if row.get("ingest_status") == MATCHED]
        if not paths:
            errors.append(f"PRICE zero-index spot-check missing: {price_key}")
        elif route == "bam" and not any("Pancreas_0.known.bed" in path for path in paths):
            errors.append(f"PRICE Pancreas_0 did not map to replicate 1: {paths}")
        elif route == "fastq" and not any("pancreas_mymapping_0.known.bed" in path for path in paths):
            errors.append(f"PRICE pancreas_mymapping_0 did not map to replicate 1: {paths}")

    pancreas_ids = {sample.sample_id for sample in EXPECTED_SAMPLE_METADATA[:6]}
    for sample_id in pancreas_ids:
        for native_class in ("annotated", "novel"):
            bam = [
                row.get("source_path", "")
                for row in rows_by_key.get(("iRibo", sample_id, "bam", native_class), [])
                if row.get("ingest_status") == MATCHED
            ]
            fastq = [
                row.get("source_path", "")
                for row in rows_by_key.get(("iRibo", sample_id, "fastq", native_class), [])
                if row.get("ingest_status") == MATCHED
            ]
            if bam and fastq and set(bam) & set(fastq):
                errors.append(f"iRibo bam/fastq paths overlap for {sample_id}/{native_class}: {bam} vs {fastq}")

    zero_counts = [
        row for row in rows
        if row.get("ingest_status") == MATCHED and int(row.get("raw_record_count") or 0) <= 0
    ]
    if zero_counts:
        errors.append(f"matched rows with zero raw_record_count: {[row.get('source_path') for row in zero_counts[:20]]}")
    ribotie_pancreas_csv = [
        row for row in rows
        if row.get("ingest_status") == MATCHED
        and row.get("tool") == "RiboTIE"
        and row.get("route") == "fastq"
        and row.get("sample_id") in pancreas_ids
    ]
    ribotie_zero = [row for row in ribotie_pancreas_csv if int(row.get("raw_record_count") or 0) <= 0]
    if ribotie_zero:
        errors.append(f"RiboTIE pancreas CSV recovery has zero-count rows: {[row.get('source_path') for row in ribotie_zero[:20]]}")

    diff_rows = [row for row in rows if row.get("reconciliation") and row.get("reconciliation") != "agree"]
    disposition_keys: set[tuple[str, str, str, str, str]] = set()
    if args.reconciliation_dispositions and args.reconciliation_dispositions.exists():
        for row in read_tsv(args.reconciliation_dispositions):
            disposition_keys.add((row.get("tool", ""), row.get("sample_id", ""), row.get("route", ""), row.get("native_class", ""), row.get("source_path", "")))
    undisposed = [
        row for row in diff_rows
        if (row.get("tool", ""), row.get("sample_id", ""), row.get("route", ""), row.get("native_class", ""), row.get("source_path", "")) not in disposition_keys
        and not row.get("reconciliation_disposition")
    ]
    if undisposed:
        errors.append(f"reconciliation diffs without disposition: {[row.get('source_path') for row in undisposed[:20]]}")

    summary_rows = [
        {"metric": "disk_leaves", "value": len(disk_leaves)},
        {"metric": "manifest_matched_or_excluded_paths", "value": len(manifest_leaf_set)},
        {"metric": "matched_rows", "value": sum(1 for row in rows if row.get("ingest_status") == MATCHED)},
        {"metric": "excluded_rows", "value": sum(1 for row in rows if row.get("ingest_status") == "excluded")},
        {"metric": "coverage_cells", "value": len(coverage_rows)},
        {"metric": "errors", "value": len(errors)},
    ]
    write_tsv(args.out_dir / "stage1_manifest_validation_summary.tsv", summary_rows, ["metric", "value"])

    if errors:
        raise SystemExit("Stage 1 manifest datachecks failed:\n  - " + "\n  - ".join(errors))
    print(f"Stage 1 manifest datachecks passed. Report: {args.out_dir}")


if __name__ == "__main__":
    main()
