#!/usr/bin/env python3
"""Stage 2 datachecks: validate translon DB fidelity against a frozen manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import sqlite3
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def manifest_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        row.get("tool", ""),
        row.get("sample_id", ""),
        row.get("input_route") or ("fastq_to_orf" if row.get("route") == "fastq" else "bam_to_orf"),
        row.get("native_class", ""),
    )


def db_rows(con: sqlite3.Connection, sql: str, params: tuple = ()) -> list[dict[str, object]]:
    con.row_factory = sqlite3.Row
    return [dict(row) for row in con.execute(sql, params)]


def table_columns(con: sqlite3.Connection, table: str) -> set[str]:
    return {row["name"] for row in db_rows(con, f"PRAGMA table_info({table})")}


def parse_bed12_records(path: Path) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    with path.open() as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            chrom, start, end, name, _score, strand = fields[:6]
            sizes = ",".join(part for part in fields[10].rstrip(",").split(",") if part)
            starts = ",".join(part for part in fields[11].rstrip(",").split(",") if part)
            records[name or f"{path.stem}_{idx}"] = {
                "bed_chrom": chrom,
                "bed_start": int(start),
                "bed_end": int(end),
                "bed_strand": strand,
                "block_sizes": sizes,
                "block_starts": starts,
            }
    return records


def parse_csv_records(path: Path) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    with path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        reader = csv.DictReader(handle, dialect=dialect)
        for idx, row in enumerate(reader, start=1):
            lower = {key.lower(): key for key in row}
            chrom_key = lower.get("chrom") or lower.get("chromosome") or lower.get("seq_region_name")
            start_key = lower.get("start") or lower.get("chrom_start") or lower.get("seq_region_start")
            end_key = lower.get("end") or lower.get("chrom_end") or lower.get("seq_region_end")
            strand_key = lower.get("strand") or lower.get("seq_region_strand")
            if not (chrom_key and start_key and end_key):
                continue
            start = int(float(row[start_key]))
            end = int(float(row[end_key]))
            if "seq_region" in start_key and start > 0:
                start -= 1
            strand = row[strand_key] if strand_key else "+"
            if strand == "1":
                strand = "+"
            elif strand == "-1":
                strand = "-"
            chrom = row[chrom_key]
            if not chrom.startswith("chr"):
                chrom = "chrM" if chrom in {"M", "MT"} else f"chr{chrom}"
            name_key = lower.get("id") or lower.get("orf_id")
            name = row.get(name_key, "") if name_key else ""
            records[name or f"{path.stem}_{idx}"] = {
                "bed_chrom": chrom,
                "bed_start": start,
                "bed_end": end,
                "bed_strand": strand,
                "block_sizes": str(end - start),
                "block_starts": "0",
            }
    return records


def raw_records(path: Path, parser: str) -> dict[str, dict[str, object]]:
    if parser == "bed12" or path.suffix == ".bed":
        return parse_bed12_records(path)
    if parser == "translonscorer_csv" or path.suffix == ".csv":
        return parse_csv_records(path)
    return {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-roundtrip-checks", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = [row for row in read_tsv(args.manifest) if row.get("ingest_status") == "matched"]
    manifest_by_path = {str(Path(row.get("source_path") or row.get("path", ""))): row for row in manifest}
    excluded_paths = {
        str(Path(row.get("source_path") or row.get("path", "")))
        for row in read_tsv(args.manifest)
        if row.get("ingest_status") == "excluded" and (row.get("source_path") or row.get("path"))
    }
    manifest_keys = {manifest_key(row) for row in manifest}

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    errors: list[str] = []
    required_parser_cols = {"path", "tool_hint", "sample_id", "input_route", "native_class", "raw_record_count", "parser_collapse_count", "translons"}
    parser_cols = table_columns(con, "parser_manifest")
    missing_cols = required_parser_cols - parser_cols
    if missing_cols:
        errors.append(f"parser_manifest missing required columns: {sorted(missing_cols)}")
    translon_cols = table_columns(con, "translons")
    if "native_class" not in translon_cols:
        errors.append("translons table missing native_class")
    if errors:
        raise SystemExit("Stage 2 DB/manifest datachecks failed:\n  - " + "\n  - ".join(errors))

    parser_rows = db_rows(con, "SELECT * FROM parser_manifest")
    parser_by_path = {str(Path(str(row["path"]))): row for row in parser_rows}
    if len(parser_by_path) != len(parser_rows):
        dupes = [path for path, count in Counter(str(Path(str(row["path"]))) for row in parser_rows).items() if count > 1]
        errors.append(f"duplicate parser_manifest paths: {dupes[:20]}")
    missing_parser = sorted(set(manifest_by_path) - set(parser_by_path))
    extra_parser = sorted(set(parser_by_path) - set(manifest_by_path))
    if missing_parser:
        errors.append(f"matched manifest paths missing in parser_manifest: {missing_parser[:20]}")
    if extra_parser:
        errors.append(f"parser_manifest paths not matched in manifest: {extra_parser[:20]}")

    count_rows: list[dict[str, object]] = []
    for path, mrow in manifest_by_path.items():
        prow = parser_by_path.get(path)
        if not prow:
            continue
        if prow["parser_status"] != "ok":
            errors.append(f"matched manifest path parsed with status {prow['parser_status']}: {path}")
        raw_count = int(mrow.get("raw_record_count") or 0)
        db_raw_count = int(prow["raw_record_count"] or 0)
        parsed = int(prow["translons"] or 0)
        collapse = int(prow["parser_collapse_count"] or 0)
        expected_parsed = db_raw_count - collapse
        count_rows.append(
            {
                "path": path,
                "raw_record_count": raw_count,
                "parser_raw_record_count": db_raw_count,
                "parser_collapse_count": collapse,
                "parser_translons": parsed,
                "status": "ok" if raw_count == db_raw_count and parsed == expected_parsed else "fail",
            }
        )
        if raw_count != db_raw_count:
            errors.append(f"raw_record_count mismatch for {path}: manifest={raw_count} parser={db_raw_count}")
        if parsed != expected_parsed:
            errors.append(f"undeclared parser collapse for {path}: raw={db_raw_count} collapse={collapse} translons={parsed}")
        if raw_count > 0 and parsed <= 0:
            errors.append(f"matched manifest path produced zero DB translons: {path}")
    write_tsv(args.out_dir / "stage2_count_conservation.tsv", count_rows, list(count_rows[0]) if count_rows else ["path"])

    translons = db_rows(con, "SELECT translon_id, source_tool, sample_id, input_route, native_class, raw_file, native_translon_id, parser_name, bed_chrom, bed_start, bed_end, bed_strand, block_count, block_sizes, block_starts, seq_region_start, seq_region_end FROM translons")
    orphan_keys = []
    orphan_paths = []
    excluded_leaks = []
    for row in translons:
        k = (str(row["source_tool"]), str(row["sample_id"]), str(row["input_route"]), str(row["native_class"]))
        raw_file = str(Path(str(row["raw_file"])))
        if k not in manifest_keys:
            orphan_keys.append({"translon_id": row["translon_id"], "key": k, "raw_file": raw_file})
        if raw_file not in manifest_by_path:
            orphan_paths.append({"translon_id": row["translon_id"], "raw_file": raw_file})
        if raw_file in excluded_paths:
            excluded_leaks.append({"translon_id": row["translon_id"], "raw_file": raw_file})
    if orphan_keys:
        errors.append(f"translons with key absent from matched manifest: {orphan_keys[:20]}")
    if orphan_paths:
        errors.append(f"translons with raw_file absent from matched manifest paths: {orphan_paths[:20]}")
    if excluded_leaks:
        errors.append(f"translons sourced from excluded manifest paths: {excluded_leaks[:20]}")

    parser_counts = {str(Path(str(row["path"]))): int(row["translons"] or 0) for row in parser_rows}
    translon_counts = Counter(str(Path(str(row["raw_file"]))) for row in translons)
    count_mismatches = [
        {"path": path, "parser_translons": parser_count, "translons": translon_counts.get(path, 0)}
        for path, parser_count in parser_counts.items()
        if parser_count != translon_counts.get(path, 0)
    ]
    if count_mismatches:
        errors.append(f"per-file parser_manifest/translons count mismatch: {count_mismatches[:20]}")

    convention_bad = db_rows(
        con,
        """
        SELECT 'translons' AS table_name, translon_id, bed_start, seq_region_start, bed_end, seq_region_end
        FROM translons
        WHERE seq_region_start != bed_start + 1 OR seq_region_end != bed_end
        UNION ALL
        SELECT 'translon_blocks' AS table_name, translon_id, bed_start, seq_region_start, bed_end, seq_region_end
        FROM translon_blocks
        WHERE seq_region_start != bed_start + 1 OR seq_region_end != bed_end
        LIMIT 20
        """,
    )
    if convention_bad:
        errors.append(f"coordinate convention failures: {convention_bad}")
    block_rows = db_rows(con, "SELECT translon_id, seq_region_strand, genomic_block_rank, translation_block_rank FROM translon_blocks ORDER BY translon_id, genomic_block_rank")
    by_tid: dict[str, list[dict[str, object]]] = {}
    for row in block_rows:
        by_tid.setdefault(str(row["translon_id"]), []).append(row)
    bad_rank = []
    for tid, blocks in by_tid.items():
        ranks = [int(row["translation_block_rank"]) for row in blocks]
        strand = int(blocks[0]["seq_region_strand"])
        if strand == 1 and ranks != sorted(ranks):
            bad_rank.append({"translon_id": tid, "strand": strand, "translation_ranks": ranks})
        if strand == -1 and ranks != sorted(ranks, reverse=True):
            bad_rank.append({"translon_id": tid, "strand": strand, "translation_ranks": ranks})
    if bad_rank:
        errors.append(f"block rank convention failures: {bad_rank[:20]}")

    roundtrip_rows: list[dict[str, object]] = []
    candidates = sorted(
        translons,
        key=lambda row: (
            0 if int(row["block_count"]) > 1 else 1,
            0 if row["bed_strand"] == "-" else 1,
            str(row["source_tool"]),
            str(row["raw_file"]),
        ),
    )
    checked = 0
    for row in candidates:
        path = Path(str(row["raw_file"]))
        manifest_row = manifest_by_path.get(str(path), {})
        records = raw_records(path, str(row["parser_name"]))
        if not records:
            continue
        raw = records.get(str(row["native_translon_id"]))
        if not raw and len(records) == 1:
            raw = next(iter(records.values()))
        if not raw:
            continue
        comparable = {
            "bed_chrom": row["bed_chrom"],
            "bed_start": int(row["bed_start"]),
            "bed_end": int(row["bed_end"]),
            "bed_strand": row["bed_strand"],
            "block_sizes": row["block_sizes"],
            "block_starts": row["block_starts"],
        }
        ok = raw == comparable
        roundtrip_rows.append(
            {
                "translon_id": row["translon_id"],
                "raw_file": str(path),
                "parser_name": row["parser_name"],
                "tool": row["source_tool"],
                "sample_id": row["sample_id"],
                "native_class": manifest_row.get("native_class", ""),
                "status": "ok" if ok else "fail",
                "raw": raw,
                "db": comparable,
            }
        )
        checked += 1
        if not ok:
            errors.append(f"raw round-trip mismatch for {row['translon_id']} from {path}: raw={raw} db={comparable}")
        if checked >= max(args.min_roundtrip_checks, 1):
            break
    write_tsv(
        args.out_dir / "stage2_raw_roundtrip.tsv",
        roundtrip_rows,
        ["translon_id", "raw_file", "parser_name", "tool", "sample_id", "native_class", "status", "raw", "db"],
    )
    if checked < args.min_roundtrip_checks:
        errors.append(f"only {checked} independent raw round-trip checks completed; required {args.min_roundtrip_checks}")

    summary_rows = [
        {"metric": "matched_manifest_rows", "value": len(manifest_by_path)},
        {"metric": "parser_manifest_rows", "value": len(parser_rows)},
        {"metric": "translons", "value": len(translons)},
        {"metric": "roundtrip_checks", "value": checked},
        {"metric": "errors", "value": len(errors)},
    ]
    write_tsv(args.out_dir / "stage2_db_manifest_validation_summary.tsv", summary_rows, ["metric", "value"])

    if errors:
        raise SystemExit("Stage 2 DB/manifest datachecks failed:\n  - " + "\n  - ".join(errors))
    print(f"Stage 2 DB/manifest datachecks passed. Report: {args.out_dir}")


if __name__ == "__main__":
    main()
