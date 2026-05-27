#!/usr/bin/env python3
"""Build sparse per-sample genomic profiles from a global unique-read BAM."""

import argparse
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

try:
    import polars as pl
except ImportError:
    pl = None

try:
    import pysam
except ImportError:
    pysam = None


READ_NAME_RE = re.compile(r"(?:^|[^\w])read_(\d+)(?:$|[^\d])")


def require_dependencies():
    if pl is None:
        raise RuntimeError("polars is required")
    if pysam is None:
        raise RuntimeError("pysam is required")


def load_manifest(path: Path) -> Tuple[Path, Dict]:
    path = Path(path)
    if path.is_dir():
        manifest_path = path / "global_matrix_manifest.json"
    else:
        manifest_path = path
    with open(manifest_path, "r") as handle:
        manifest = json.load(handle)
    if manifest.get("matrix_format") != "sparse-parquet":
        raise ValueError(f"{manifest_path} is not a sparse-parquet manifest")
    return manifest_path, manifest


def generation_paths(manifest_path: Path, manifest: Dict) -> List[Path]:
    generations = manifest.get("generations") or [
        {
            "generation": manifest.get("generation", 1),
            "path": manifest.get("path", "global_counts/generation=000001"),
        }
    ]
    return [manifest_path.parent / item["path"] for item in generations]


def parse_read_id(query_name: str) -> Optional[int]:
    match = READ_NAME_RE.search(query_name)
    if match is None:
        return None
    return int(match.group(1))


def alignment_event(record) -> Optional[Tuple[str, int, str]]:
    if record.is_unmapped:
        return None
    chrom = record.reference_name
    if chrom is None:
        return None
    strand = "-" if record.is_reverse else "+"
    return chrom, int(record.reference_start), strand


def load_tombstone_filters(manifest_path: Path, manifest: Dict) -> Tuple[Set[int], Set[str], Set[str]]:
    tombstones_path = manifest_path.parent / manifest.get("tombstones_path", "global_tombstones.parquet")
    if not tombstones_path.exists():
        return set(), set(), set()
    frame = pl.read_parquet(tombstones_path)
    if frame.is_empty():
        return set(), set(), set()
    frame = frame.filter(pl.col("active") == True)
    read_ids = set()
    sample_names = set()
    study_ids = set()
    for row in frame.to_dicts():
        if row.get("read_id") is not None:
            read_ids.add(int(row["read_id"]))
        if row.get("sample_id") is not None:
            sample_names.add(str(row["sample_id"]))
        if row.get("study_id") is not None:
            study_ids.add(str(row["study_id"]))
    return read_ids, sample_names, study_ids


class BucketEventSpool:
    """Spool BAM alignment events to bounded per-read-bucket TSV files."""

    def __init__(self, temp_dir: Path):
        self.temp_dir = Path(temp_dir)
        self.handles = {}
        self.paths = {}
        self.rows = 0
        self.skipped_unparsed = 0
        self.skipped_unmapped = 0

    def write(self, read_bucket: int, read_id: int, chrom: str, position: int, strand: str):
        handle = self.handles.get(read_bucket)
        if handle is None:
            path = self.temp_dir / f"events.read_bucket={read_bucket:06d}.tsv"
            handle = open(path, "w")
            handle.write("read_bucket\tread_id\tchrom\tposition\tstrand\tevent_count\n")
            self.handles[read_bucket] = handle
            self.paths[read_bucket] = path
        handle.write(f"{read_bucket}\t{read_id}\t{chrom}\t{position}\t{strand}\t1\n")
        self.rows += 1

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def spool_bam_events(
    bam_path: Path,
    read_bucket_size: int,
    temp_dir: Path,
    exclude_read_ids: Optional[Set[int]] = None,
) -> BucketEventSpool:
    spool = BucketEventSpool(temp_dir)
    exclude_read_ids = exclude_read_ids or set()
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for record in bam.fetch(until_eof=True):
                read_id = parse_read_id(record.query_name)
                if read_id is None:
                    spool.skipped_unparsed += 1
                    continue
                if read_id in exclude_read_ids:
                    continue
                event = alignment_event(record)
                if event is None:
                    spool.skipped_unmapped += 1
                    continue
                chrom, position, strand = event
                spool.write(read_id // read_bucket_size, read_id, chrom, position, strand)
    finally:
        spool.close()
    return spool


def scan_event_bucket(path: Path):
    schema = {
        "read_bucket": pl.UInt32,
        "read_id": pl.UInt64,
        "chrom": pl.Utf8,
        "position": pl.UInt32,
        "strand": pl.Utf8,
        "event_count": pl.UInt32,
    }
    try:
        return pl.scan_csv(path, separator="\t", schema_overrides=schema)
    except TypeError:
        return pl.scan_csv(path, separator="\t", dtypes=schema)


def count_bucket_paths(generation_roots: Iterable[Path], read_bucket: int) -> List[Path]:
    paths = []
    bucket_name = f"read_bucket={read_bucket:06d}"
    for root in generation_roots:
        paths.extend(sorted((root / bucket_name).glob("*.parquet")))
    return paths


def load_count_bucket(generation_roots: Iterable[Path], read_bucket: int):
    paths = count_bucket_paths(generation_roots, read_bucket)
    if not paths:
        return None
    return pl.scan_parquet([str(path) for path in paths])


def load_samples(manifest_path: Path, manifest: Dict):
    samples_path = manifest_path.parent / manifest.get("lookup_tables", {}).get("samples", "global_samples.parquet")
    if not samples_path.exists():
        raise FileNotFoundError(f"Sample lookup table does not exist: {samples_path}")
    return pl.scan_parquet(str(samples_path)).select(
        "sample_id",
        "sample_name",
        "study_id",
        "study_id_int",
    )


def filter_count_facts(counts, samples, sample_tombstones: Set[str], study_tombstones: Set[str]):
    if not sample_tombstones and not study_tombstones:
        return counts
    annotated = counts.join(samples, on=["sample_id", "study_id_int"], how="left")
    if sample_tombstones:
        annotated = annotated.filter(~pl.col("sample_name").is_in(sorted(sample_tombstones)))
    if study_tombstones:
        annotated = annotated.filter(~pl.col("study_id").is_in(sorted(study_tombstones)))
    return annotated.select("read_bucket", "read_id", "sample_id", "study_id_int", "count")


def write_bucket_profile(
    event_path: Path,
    read_bucket: int,
    generation_roots: Iterable[Path],
    samples,
    sample_tombstones: Set[str],
    study_tombstones: Set[str],
    output_path: Path,
) -> int:
    counts = load_count_bucket(generation_roots, read_bucket)
    if counts is None:
        return 0

    events = scan_event_bucket(event_path)
    counts = filter_count_facts(counts, samples, sample_tombstones, study_tombstones)
    profile = (
        events.join(counts, on="read_id", how="inner")
        .with_columns((pl.col("event_count").cast(pl.UInt64) * pl.col("count").cast(pl.UInt64)).alias("weighted_count"))
        .group_by("sample_id", "study_id_int", "chrom", "position", "strand")
        .agg(pl.col("weighted_count").sum().alias("count"))
        .select(
            pl.col("sample_id").cast(pl.UInt32),
            pl.col("study_id_int").cast(pl.UInt32),
            "chrom",
            pl.col("position").cast(pl.UInt32),
            "strand",
            pl.col("count").cast(pl.UInt64),
        )
    )
    frame = profile.collect()
    if frame.is_empty():
        return 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.write_parquet(output_path, compression="zstd")
    return frame.height


def coalesce_profile_shards(shards_dir: Path, output: Path):
    paths = sorted(shards_dir.glob("*.parquet"))
    if not paths:
        pl.DataFrame(
            {
                "sample_id": pl.Series([], dtype=pl.UInt32),
                "study_id_int": pl.Series([], dtype=pl.UInt32),
                "chrom": pl.Series([], dtype=pl.Utf8),
                "position": pl.Series([], dtype=pl.UInt32),
                "strand": pl.Series([], dtype=pl.Utf8),
                "count": pl.Series([], dtype=pl.UInt64),
            }
        ).write_parquet(output)
        return 0

    output.parent.mkdir(parents=True, exist_ok=True)
    scan = pl.scan_parquet([str(path) for path in paths])
    result = (
        scan.group_by("sample_id", "study_id_int", "chrom", "position", "strand")
        .agg(pl.col("count").sum().alias("count"))
        .sort("sample_id", "chrom", "position", "strand")
        .collect()
    )
    result.write_parquet(output, compression="zstd")
    return result.height


def build_profile(
    bam: Path,
    manifest: Path,
    output: Path,
    keep_bucket_shards: bool,
    no_final_coalesce: bool,
):
    require_dependencies()
    manifest_path, manifest_data = load_manifest(manifest)
    read_bucket_size = int(manifest_data["read_bucket_size"])
    generation_roots = generation_paths(manifest_path, manifest_data)
    samples = load_samples(manifest_path, manifest_data)
    read_tombstones, sample_tombstones, study_tombstones = load_tombstone_filters(manifest_path, manifest_data)

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="profile_global_bam_", dir=output.parent) as temp_name:
        temp_dir = Path(temp_name)
        events_dir = temp_dir / "events"
        shards_dir = temp_dir / "profile_shards"
        events_dir.mkdir()
        shards_dir.mkdir()

        event_spool = spool_bam_events(bam, read_bucket_size, events_dir, read_tombstones)
        rows_written = 0
        buckets_loaded = 0
        for read_bucket, event_path in sorted(event_spool.paths.items()):
            count_paths = count_bucket_paths(generation_roots, read_bucket)
            if not count_paths:
                continue
            buckets_loaded += 1
            shard_path = shards_dir / f"profile.read_bucket={read_bucket:06d}.parquet"
            rows_written += write_bucket_profile(
                event_path,
                read_bucket,
                generation_roots,
                samples,
                sample_tombstones,
                study_tombstones,
                shard_path,
            )

        if no_final_coalesce:
            if output.exists():
                if output.is_dir():
                    shutil.rmtree(output)
                else:
                    output.unlink()
            shutil.copytree(shards_dir, output)
            final_rows = rows_written
        else:
            final_rows = coalesce_profile_shards(shards_dir, output)
            if keep_bucket_shards:
                bucket_output = output.with_suffix(output.suffix + ".buckets")
                if bucket_output.exists():
                    shutil.rmtree(bucket_output)
                shutil.copytree(shards_dir, bucket_output)

        return {
            "bam_events": event_spool.rows,
            "skipped_unparsed": event_spool.skipped_unparsed,
            "skipped_unmapped": event_spool.skipped_unmapped,
            "touched_buckets": len(event_spool.paths),
            "loaded_count_buckets": buckets_loaded,
            "bucket_profile_rows": rows_written,
            "final_profile_rows": final_rows,
        }


def main():
    parser = argparse.ArgumentParser(
        description="Build sparse per-sample genomic profiles from a global unique-read BAM"
    )
    parser.add_argument("--bam", required=True, type=Path, help="Global unique-read BAM")
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Sparse global matrix manifest or output directory",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output sparse profile Parquet file, or directory with --no-final-coalesce",
    )
    parser.add_argument(
        "--keep-bucket-shards",
        action="store_true",
        help="Keep intermediate per-read-bucket profile shards beside the final output",
    )
    parser.add_argument(
        "--no-final-coalesce",
        action="store_true",
        help="Write per-bucket profile shards only; downstream consumers must sum duplicate loci across buckets",
    )
    args = parser.parse_args()

    try:
        stats = build_profile(
            args.bam,
            args.manifest,
            args.output,
            args.keep_bucket_shards,
            args.no_final_coalesce,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    print(json.dumps(stats, sort_keys=True), file=sys.stderr)


if __name__ == "__main__":
    main()
