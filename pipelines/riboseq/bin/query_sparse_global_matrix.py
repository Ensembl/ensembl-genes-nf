#!/usr/bin/env python3
"""Query and maintain sparse global unique-read count stores."""

import argparse
import json
import os
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

try:
    import polars as pl
except ImportError:
    pl = None


def require_polars():
    if pl is None:
        raise RuntimeError("polars is required to query sparse global matrix stores")


def parse_csv_ints(value: Optional[str]) -> Optional[List[int]]:
    if not value:
        return None
    return [int(item) for item in value.split(",") if item != ""]


def parse_csv_strings(value: Optional[str]) -> Optional[List[str]]:
    if not value:
        return None
    return [item for item in value.split(",") if item != ""]


def load_manifest(manifest_path: Path) -> Dict:
    with open(manifest_path, "r") as handle:
        manifest = json.load(handle)
    if manifest.get("matrix_format") != "sparse-parquet":
        raise ValueError(f"{manifest_path} is not a sparse-parquet global matrix manifest")
    return manifest


def resolve_paths(path: Path) -> tuple[Path, Path, Dict]:
    path = Path(path)
    if path.is_dir() and (path / "global_matrix_manifest.json").exists():
        manifest_path = path / "global_matrix_manifest.json"
    elif path.name.endswith("_matrix_manifest.json") or path.name == "global_matrix_manifest.json":
        manifest_path = path
    else:
        raise ValueError(
            "Expected an output directory containing global_matrix_manifest.json "
            "or a direct manifest path"
        )

    manifest = load_manifest(manifest_path)
    matrix_path = manifest_path.parent / manifest["path"]
    if not matrix_path.exists():
        raise FileNotFoundError(f"Sparse matrix path from manifest does not exist: {matrix_path}")
    return manifest_path, matrix_path, manifest


def generation_paths(manifest_path: Path, manifest: Dict) -> List[Path]:
    generations = manifest.get("generations") or [
        {
            "generation": manifest.get("generation", 1),
            "path": manifest.get("path", "global_matrix.parquet"),
        }
    ]
    return [manifest_path.parent / item["path"] for item in generations]


def empty_query_frame():
    return pl.DataFrame(
        {
            "read_id": pl.Series([], dtype=pl.UInt64),
            "sample_index": pl.Series([], dtype=pl.UInt32),
            "sample_id": pl.Series([], dtype=pl.Utf8),
            "study_id": pl.Series([], dtype=pl.Utf8),
            "count": pl.Series([], dtype=pl.UInt32),
        }
    )


def project_query_columns(scan, manifest_path: Path, manifest: Dict):
    schema = manifest.get("schema", {})
    if schema.get("sample_id") == "uint32":
        lookup_tables = manifest.get("lookup_tables", {})
        samples_path = manifest_path.parent / lookup_tables.get("samples", "global_samples.parquet")
        if not samples_path.exists():
            raise FileNotFoundError(f"Sparse sample lookup table does not exist: {samples_path}")
        samples = pl.scan_parquet(str(samples_path)).select(
            "sample_id",
            "sample_name",
            "study_id",
        )
        return (
            scan.join(samples, on="sample_id", how="left")
            .select(
                "read_id",
                pl.col("sample_id").alias("sample_index"),
                pl.col("sample_name").alias("sample_id"),
                "study_id",
                "count",
            )
        )

    return scan.select("read_id", "sample_index", "sample_id", "study_id", "count")


def tombstone_filter_expr(tombstones: List[Dict]):
    expr = None
    for tombstone in tombstones:
        if tombstone.get("active") is False:
            continue

        current = None
        if tombstone.get("read_id") is not None:
            current = pl.col("read_id") == int(tombstone["read_id"])
        if tombstone.get("sample_id") is not None:
            sample_expr = pl.col("sample_id") == str(tombstone["sample_id"])
            current = sample_expr if current is None else current & sample_expr
        if tombstone.get("study_id") is not None:
            study_expr = pl.col("study_id") == str(tombstone["study_id"])
            current = study_expr if current is None else current & study_expr

        if current is not None:
            expr = current if expr is None else expr | current

    return expr


def load_tombstones(manifest_path: Path, manifest: Dict) -> List[Dict]:
    tombstone_path = manifest_path.parent / manifest.get("tombstones_path", "global_tombstones.parquet")
    legacy = manifest.get("tombstones", [])
    if not tombstone_path.exists():
        return legacy
    frame = pl.read_parquet(tombstone_path)
    if frame.is_empty():
        return legacy
    return legacy + frame.to_dicts()


def query_sparse_store(
    store: Path,
    read_ids: Optional[List[int]],
    sample_ids: Optional[List[str]],
    study_ids: Optional[List[str]],
    include_tombstoned: bool,
):
    require_polars()
    manifest_path, matrix_path, manifest = resolve_paths(store)
    matrix_paths = generation_paths(manifest_path, manifest)

    if read_ids is not None:
        bucket_size = int(manifest["read_bucket_size"])
        bucket_paths = []
        for matrix_path_item in matrix_paths:
            for bucket in sorted({read_id // bucket_size for read_id in read_ids}):
                bucket_paths.extend(sorted((matrix_path_item / f"read_bucket={bucket:06d}").glob("*.parquet")))
        if not bucket_paths:
            return empty_query_frame()
        scan = pl.scan_parquet([str(path) for path in bucket_paths])
    else:
        scan = pl.concat([
            pl.scan_parquet(str(path / "read_bucket=*" / "*.parquet"))
            for path in matrix_paths
        ])

    scan = project_query_columns(scan, manifest_path, manifest)

    if read_ids is not None:
        scan = scan.filter(pl.col("read_id").is_in(read_ids))
    if sample_ids is not None:
        scan = scan.filter(pl.col("sample_id").is_in(sample_ids))
    if study_ids is not None:
        scan = scan.filter(pl.col("study_id").is_in(study_ids))

    if not include_tombstoned:
        tombstone_expr = tombstone_filter_expr(load_tombstones(manifest_path, manifest))
        if tombstone_expr is not None:
            scan = scan.filter(~tombstone_expr)

    return (
        scan.select("read_id", "sample_index", "sample_id", "study_id", "count")
        .sort("read_id", "sample_index", "study_id")
        .collect()
    )


def atomic_write_json(path: Path, payload: Dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        shutil.move(str(tmp_path), path)
    finally:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass


def atomic_write_parquet(path: Path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        frame.write_parquet(tmp_path)
        shutil.move(str(tmp_path), path)
    finally:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass


def add_tombstone(
    store: Path,
    reason: str,
    read_id: Optional[int],
    sample_id: Optional[str],
    study_id: Optional[str],
):
    manifest_path, _, manifest = resolve_paths(store)
    if read_id is None and sample_id is None and study_id is None:
        raise ValueError("At least one tombstone selector is required")

    tombstone_path = manifest_path.parent / manifest.get("tombstones_path", "global_tombstones.parquet")
    if tombstone_path.exists():
        current = pl.read_parquet(tombstone_path)
    else:
        current = pl.DataFrame(
            {
                "tombstone_id": pl.Series([], dtype=pl.UInt64),
                "active": pl.Series([], dtype=pl.Boolean),
                "created_at": pl.Series([], dtype=pl.Utf8),
                "reason": pl.Series([], dtype=pl.Utf8),
                "read_id": pl.Series([], dtype=pl.UInt64),
                "sample_id": pl.Series([], dtype=pl.Utf8),
                "study_id": pl.Series([], dtype=pl.Utf8),
            }
        )

    next_id = 0 if current.is_empty() else int(current["tombstone_id"].max()) + 1
    row = pl.DataFrame(
        {
            "tombstone_id": pl.Series([next_id], dtype=pl.UInt64),
            "active": pl.Series([True], dtype=pl.Boolean),
            "created_at": pl.Series([datetime.now(timezone.utc).isoformat()], dtype=pl.Utf8),
            "reason": pl.Series([reason], dtype=pl.Utf8),
            "read_id": pl.Series([read_id], dtype=pl.UInt64),
            "sample_id": pl.Series([sample_id], dtype=pl.Utf8),
            "study_id": pl.Series([study_id], dtype=pl.Utf8),
        }
    )
    atomic_write_parquet(tombstone_path, pl.concat([current, row]))
    manifest["tombstones_path"] = str(tombstone_path.relative_to(manifest_path.parent))
    manifest.pop("tombstones", None)
    atomic_write_json(manifest_path, manifest)
    return tombstone_path


def materialize_retained_view(store: Path, output: Optional[Path]):
    frame = query_sparse_store(
        store,
        read_ids=None,
        sample_ids=None,
        study_ids=None,
        include_tombstoned=False,
    )
    manifest_path, _, manifest = resolve_paths(store)
    if output is None:
        output = manifest_path.parent / manifest.get("retained_counts", "global_retained_counts.parquet")
    raw = frame.select(
        "read_id",
        pl.col("sample_index").alias("sample_id"),
        "study_id",
        "count",
    )
    atomic_write_parquet(output, raw)
    return output


def main():
    parser = argparse.ArgumentParser(
        description="Query and maintain sparse global unique-read matrix stores"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    query = subparsers.add_parser("query", help="Query sparse count facts")
    query.add_argument("store", type=Path, help="Global output directory or manifest path")
    query.add_argument("--read-ids", help="Comma-separated global read IDs")
    query.add_argument("--sample-ids", help="Comma-separated sample IDs")
    query.add_argument("--study-ids", help="Comma-separated study IDs")
    query.add_argument("--include-tombstoned", action="store_true")
    query.add_argument("--output", type=Path, help="Write TSV to this path instead of stdout")

    tombstone = subparsers.add_parser("tombstone", help="Add a logical tombstone to the manifest")
    tombstone.add_argument("store", type=Path, help="Global output directory or manifest path")
    tombstone.add_argument("--read-id", type=int)
    tombstone.add_argument("--sample-id")
    tombstone.add_argument("--study-id")
    tombstone.add_argument("--reason", required=True)

    retained = subparsers.add_parser("retained-view", help="Materialize a retained-count Parquet view")
    retained.add_argument("store", type=Path, help="Global output directory or manifest path")
    retained.add_argument("--output", type=Path, help="Output parquet path")

    args = parser.parse_args()

    try:
        if args.command == "query":
            frame = query_sparse_store(
                args.store,
                parse_csv_ints(args.read_ids),
                parse_csv_strings(args.sample_ids),
                parse_csv_strings(args.study_ids),
                args.include_tombstoned,
            )
            if args.output:
                frame.write_csv(args.output, separator="\t")
            else:
                sys.stdout.write(frame.write_csv(separator="\t"))
        elif args.command == "tombstone":
            path = add_tombstone(
                args.store,
                args.reason,
                args.read_id,
                args.sample_id,
                args.study_id,
            )
            print(f"Updated tombstones in {path}", file=sys.stderr)
        elif args.command == "retained-view":
            path = materialize_retained_view(args.store, args.output)
            print(f"Wrote retained view to {path}", file=sys.stderr)
        else:
            raise AssertionError(args.command)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
