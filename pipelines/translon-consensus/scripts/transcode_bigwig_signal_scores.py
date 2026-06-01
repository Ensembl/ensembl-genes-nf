#!/usr/bin/env python3
"""Score translon DB features from stranded BigWig signal.

The scorer is deliberately sample-oriented: each sample is loaded, deduplicated,
scored, written to a checkpoint file, and then released.  This keeps peak memory
bounded while still avoiding repeated BigWig extraction for features called by
multiple tools.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
import sqlite3
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    import pyBigWig

    HAS_PYBIGWIG = True
except ImportError:  # pragma: no cover - exercised in notebook environments
    pyBigWig = None
    HAS_PYBIGWIG = False


KEY_COLS = [
    "feature_key",
    "sample_id",
    "bed_chrom",
    "bed_start",
    "bed_end",
    "bed_strand",
    "block_sizes",
    "block_starts",
    "spliced_length_nt",
]


@dataclass(frozen=True)
class ScoreConfig:
    """Runtime and metric controls.

    `max_features_per_sample` is the fastest way to get usable signal across all
    samples first.  Set it to None for full production scoring once the pilot
    output looks sane.
    """

    flank_nt: int = 30
    body_edge_nt: int = 15
    min_body_nt: int = 30
    min_total_signal: float = 10.0
    max_features_per_sample: int | None = 2000
    max_total_features: int | None = None
    psite_offset: int = 0
    overwrite: bool = False


def discover_bigwigs(root: Path) -> pd.DataFrame:
    """Discover `{sample}.{forward,reverse}.bw` style stranded BigWigs."""

    if not root.exists():
        return pd.DataFrame(columns=["sample_id", "fwd_path", "rev_path", "has_pair"])

    suffixes = {
        "fwd_path": (".forward.bw", ".fwd.bw", "_forward.bw"),
        "rev_path": (".reverse.bw", ".rev.bw", "_reverse.bw"),
    }
    maps: dict[str, dict[str, Path]] = {"fwd_path": {}, "rev_path": {}}

    for bw in root.rglob("*.bw"):
        name = bw.name
        for col, col_suffixes in suffixes.items():
            for suffix in col_suffixes:
                if name.endswith(suffix):
                    maps[col][name[: -len(suffix)]] = bw
                    break
            else:
                continue
            break

    rows = []
    for sample_id in sorted(set(maps["fwd_path"]) | set(maps["rev_path"])):
        fwd = maps["fwd_path"].get(sample_id)
        rev = maps["rev_path"].get(sample_id)
        rows.append(
            {
                "sample_id": sample_id,
                "fwd_path": str(fwd) if fwd else None,
                "rev_path": str(rev) if rev else None,
                "has_pair": fwd is not None and rev is not None,
            }
        )
    return pd.DataFrame(rows)


def parse_bed12_blocks(
    bed_start: int, block_sizes_str: str, block_starts_str: str
) -> list[tuple[int, int]]:
    sizes = [int(x) for x in str(block_sizes_str).rstrip(",").split(",") if x]
    starts = [int(x) for x in str(block_starts_str).rstrip(",").split(",") if x]
    if len(sizes) != len(starts):
        raise ValueError("block_sizes and block_starts have different lengths")
    intervals = [(bed_start + start, bed_start + start + size) for start, size in zip(starts, sizes)]
    return sorted(intervals)


def _safe_ratio(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator):
        return np.nan
    if denominator <= 0:
        return np.nan if numerator <= 0 else np.inf
    return float(numerator / denominator)


def _mean_positive_window(arr: np.ndarray | None) -> tuple[float, int, float]:
    if arr is None or len(arr) == 0:
        return np.nan, 0, 0.0
    clean = np.asarray(arr, dtype=np.float32)
    clean = clean[np.isfinite(clean)]
    if len(clean) == 0:
        return np.nan, 0, 0.0
    return float(clean.mean()), int(np.count_nonzero(clean > 0)), float(clean.sum())


def body_slice(arr: np.ndarray, edge_nt: int, min_body_nt: int) -> tuple[np.ndarray, int, int]:
    start = min(edge_nt, len(arr))
    end = max(start, len(arr) - edge_nt)
    if end - start < min_body_nt:
        start, end = 0, len(arr)
    return arr[start:end], start, end


def mean_coverage(arr: np.ndarray | None) -> float:
    if arr is None or len(arr) == 0:
        return np.nan
    return float(np.mean(arr))


def periodicity_score(arr, start_offset=0):
    """Fraction of body signal in the annotated ORF frame."""
    if arr is None or len(arr) < 9:
        return np.nan
    positions = (np.arange(len(arr)) + start_offset) % 3
    frame_sums = np.array([arr[positions == f].sum() for f in range(3)], float)
    total = frame_sums.sum()
    if total <= 0:
        return np.nan
    return float(frame_sums[0] / total)


def uniformity_score(arr: np.ndarray | None) -> float:
    """Codon-binned body uniformity, `1 / (1 + CV)`, in the range 0..1."""

    if arr is None or len(arr) < 9:
        return np.nan
    trimmed = arr[: len(arr) - (len(arr) % 3)]
    if len(trimmed) < 9:
        return np.nan
    codon_signal = trimmed.reshape(-1, 3).sum(axis=1).astype(float)
    mean = codon_signal.mean()
    if mean <= 0:
        return np.nan
    return float(1.0 / (1.0 + codon_signal.std(ddof=0) / mean))


class BigWigPair:
    def __init__(self, fwd_path: str, rev_path: str):
        if not HAS_PYBIGWIG:
            raise RuntimeError("pyBigWig is not installed")
        self.handles = {
            "+": pyBigWig.open(str(fwd_path)),
            "-": pyBigWig.open(str(rev_path)),
        }

    def close(self) -> None:
        for handle in self.handles.values():
            try:
                handle.close()
            except Exception:
                pass

    def chrom_len(self, strand: str, chrom: str) -> int | None:
        chroms = self.handles[strand].chroms()
        return chroms.get(chrom)

    def values(self, chrom: str, start: int, end: int, strand: str, reverse: bool = False) -> np.ndarray | None:
        chrom_len = self.chrom_len(strand, chrom)
        if chrom_len is None:
            return None
        start = max(0, int(start))
        end = min(int(end), chrom_len)
        if start >= end:
            return None
        vals = self.handles[strand].values(chrom, start, end, numpy=True)
        if vals is None:
            arr = np.zeros(end - start, dtype=np.float32)
        else:
            arr = np.where(np.isnan(vals), 0.0, vals).astype(np.float32)
        return arr[::-1] if reverse else arr

    def transcript_values(self, chrom: str, intervals: Iterable[tuple[int, int]], strand: str) -> np.ndarray | None:
        chunks = []
        for start, end in intervals:
            vals = self.values(chrom, start, end, strand, reverse=False)
            if vals is not None:
                chunks.append(vals)
        if not chunks:
            return None
        arr = np.concatenate(chunks)
        return arr[::-1] if strand == "-" else arr

    def upstream_of_start(self, chrom: str, intervals: list[tuple[int, int]], strand: str, nt: int) -> np.ndarray | None:
        if strand == "+":
            boundary = intervals[0][0]
            return self.values(chrom, boundary - nt, boundary, strand)
        boundary = intervals[-1][1]
        return self.values(chrom, boundary, boundary + nt, strand, reverse=True)

    def downstream_of_stop(self, chrom: str, intervals: list[tuple[int, int]], strand: str, nt: int) -> np.ndarray | None:
        if strand == "+":
            boundary = intervals[-1][1]
            return self.values(chrom, boundary, boundary + nt, strand)
        boundary = intervals[0][0]
        return self.values(chrom, boundary - nt, boundary, strand, reverse=True)


def score_signal_array(arr: np.ndarray | None, cfg: ScoreConfig) -> dict[str, float]:
    if arr is None or len(arr) == 0:
        return {
            "mean_cov": np.nan,
            "body_mean_cov": np.nan,
            "periodicity": np.nan,
            "uniformity": np.nan,
            "body_total_signal": 0.0,
            "body_nonzero_nt": 0,
            "coverage_pass": False,
        }

    body, body_start, _body_end = body_slice(arr, cfg.body_edge_nt, cfg.min_body_nt)
    body_mean, body_nonzero, body_total = _mean_positive_window(body)
    return {
        "mean_cov": mean_coverage(arr),
        "body_mean_cov": body_mean,
        "periodicity": periodicity_score(body, start_offset=body_start + cfg.psite_offset),
        "uniformity": uniformity_score(body),
        "body_total_signal": body_total,
        "body_nonzero_nt": body_nonzero,
        "coverage_pass": bool(body_total >= cfg.min_total_signal),
    }


def score_feature(row: pd.Series, bws: BigWigPair, cfg: ScoreConfig) -> dict[str, object]:
    chrom = str(row["bed_chrom"])
    strand = str(row["bed_strand"])
    intervals = parse_bed12_blocks(int(row["bed_start"]), row["block_sizes"], row["block_starts"])
    arr = bws.transcript_values(chrom, intervals, strand)

    metrics: dict[str, object] = {
        "feature_key": row["feature_key"],
        "sample_id": row["sample_id"],
        "bed_chrom": chrom,
        "bed_start": int(row["bed_start"]),
        "bed_end": int(row["bed_end"]),
        "bed_strand": strand,
        "spliced_length_nt": int(row["spliced_length_nt"]),
        "extract_status": "ok" if arr is not None else "no_orf_signal",
        "orf_signal_nt": 0 if arr is None else int(len(arr)),
    }
    metrics.update(score_signal_array(arr, cfg))

    if arr is None:
        return metrics

    start_window = arr[: cfg.flank_nt]
    pre_stop_window = arr[max(0, len(arr) - cfg.flank_nt) :]
    upstream = bws.upstream_of_start(chrom, intervals, strand, cfg.flank_nt)
    downstream = bws.downstream_of_stop(chrom, intervals, strand, cfg.flank_nt)

    start_mean, start_nonzero, start_total = _mean_positive_window(start_window)
    upstream_mean, upstream_nonzero, upstream_total = _mean_positive_window(upstream)
    pre_stop_mean, pre_stop_nonzero, pre_stop_total = _mean_positive_window(pre_stop_window)
    post_stop_mean, post_stop_nonzero, post_stop_total = _mean_positive_window(downstream)

    metrics.update(
        {
            "start_window_mean": start_mean,
            "start_upstream_mean": upstream_mean,
            "start_rise_ratio": _safe_ratio(start_mean, upstream_mean),
            "stop_pre_mean": pre_stop_mean,
            "stop_post_mean": post_stop_mean,
            "stop_dropoff_ratio": _safe_ratio(post_stop_mean, pre_stop_mean),
            "start_window_nonzero_nt": start_nonzero,
            "start_upstream_nonzero_nt": upstream_nonzero,
            "stop_pre_nonzero_nt": pre_stop_nonzero,
            "stop_post_nonzero_nt": post_stop_nonzero,
            "start_window_total_signal": start_total,
            "start_upstream_total_signal": upstream_total,
            "stop_pre_total_signal": pre_stop_total,
            "stop_post_total_signal": post_stop_total,
            "upstream_flank_nt": 0 if upstream is None else int(len(upstream)),
            "downstream_flank_nt": 0 if downstream is None else int(len(downstream)),
            "flank_note": "genomic_flanks_not_splice_aware",
        }
    )
    return metrics


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def _read_sql(con: sqlite3.Connection, sql: str, params: tuple[object, ...] = ()) -> pd.DataFrame:
    return pd.read_sql_query(sql, con, params=params)


def load_psite_offsets(path: Path | None) -> dict[str, int]:
    """Load frozen per-sample mod-3 P-site offsets."""

    if path is None:
        logging.warning("No --psite-offsets path provided; using offset 0 for all samples")
        return {}
    if not path.exists():
        logging.warning("P-site offsets file %s is absent; using offset 0 for all samples", path)
        return {}

    raw = json.loads(path.read_text())
    offsets: dict[str, int] = {}
    for sample_id, value in raw.items():
        offset = value.get("offset") if isinstance(value, dict) else value
        try:
            offsets[str(sample_id)] = int(offset) % 3
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid P-site offset for {sample_id!r}: {offset!r}") from exc
    return offsets


def load_sample_orfs(con: sqlite3.Connection, sample_id: str, cfg: ScoreConfig) -> pd.DataFrame:
    limit = f"LIMIT {int(cfg.max_features_per_sample)}" if cfg.max_features_per_sample else ""
    return _read_sql(
        con,
        f"""
        SELECT DISTINCT
               source_tool,
               source_feature_class AS cls,
               feature_key,
               bed_chrom,
               bed_start,
               bed_end,
               bed_strand,
               block_sizes,
               block_starts,
               spliced_length_nt,
               sample_id
        FROM translons
        WHERE qc_status = 'pass'
          AND source_feature_class IN ('cds', 'non_cds')
          AND sample_id = ?
          AND sample_id NOT GLOB '*_fastq'
        ORDER BY feature_key, source_tool
        {limit}
        """,
        (sample_id,),
    )


def score_database(
    db: Path,
    bigwig_root: Path,
    out_dir: Path,
    cfg: ScoreConfig,
    psite_offsets: dict[str, int] | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    if not HAS_PYBIGWIG:
        raise RuntimeError("pyBigWig is required for signal extraction")
    if not db.exists():
        raise FileNotFoundError(db)

    out_dir.mkdir(parents=True, exist_ok=True)
    batch_dir = out_dir / "signal_score_batches"
    batch_dir.mkdir(exist_ok=True)

    bw_manifest = discover_bigwigs(bigwig_root)
    bw_manifest = bw_manifest[bw_manifest["has_pair"]].copy()
    if bw_manifest.empty:
        raise FileNotFoundError(f"No paired BigWigs found under {bigwig_root}")

    psite_offsets = psite_offsets or {}
    summaries: list[dict[str, object]] = []
    total_unique_scored = 0
    con = sqlite3.connect(db)
    try:
        available_samples = set(
            _read_sql(
                con,
                """
                SELECT DISTINCT sample_id FROM translons
                WHERE qc_status = 'pass'
                  AND source_feature_class IN ('cds', 'non_cds')
                  AND sample_id NOT GLOB '*_fastq'
                """,
            )["sample_id"]
        )
        sample_rows = bw_manifest[bw_manifest["sample_id"].isin(available_samples)].sort_values("sample_id")

        for _, bw_row in sample_rows.iterrows():
            sample_id = str(bw_row["sample_id"])
            sample_cfg = replace(cfg, psite_offset=psite_offsets.get(sample_id, 0))
            done = batch_dir / f"{_slug(sample_id)}.done.json"
            out_path = batch_dir / f"{_slug(sample_id)}.scores.tsv.gz"
            if done.exists() and out_path.exists() and not cfg.overwrite:
                previous_summary = json.loads(done.read_text())
                if previous_summary.get("psite_offset") == sample_cfg.psite_offset:
                    summaries.append(previous_summary)
                    continue
                logging.info(
                    "Sample %s: existing batch lacks matching P-site offset; recomputing",
                    sample_id,
                )
            if cfg.max_total_features is not None and total_unique_scored >= cfg.max_total_features:
                break
            logging.info("Sample %s: applying P-site offset %s", sample_id, sample_cfg.psite_offset)

            sample_orfs = load_sample_orfs(con, sample_id, sample_cfg)
            if sample_orfs.empty:
                continue
            unique_features = sample_orfs[KEY_COLS].drop_duplicates().reset_index(drop=True)
            if cfg.max_total_features is not None:
                remaining = max(0, cfg.max_total_features - total_unique_scored)
                unique_features = unique_features.head(remaining)

            annotations = (
                sample_orfs.groupby(["feature_key", "sample_id"], sort=False)
                .agg(
                    source_tools=("source_tool", lambda s: "|".join(sorted(set(map(str, s))))),
                    classes=("cls", lambda s: "|".join(sorted(set(map(str, s))))),
                    tool_rows=("source_tool", "size"),
                )
                .reset_index()
            )
            annotation_lookup = {
                (row.feature_key, row.sample_id): row
                for row in annotations.itertuples(index=False)
            }

            records = []
            bws = BigWigPair(str(bw_row["fwd_path"]), str(bw_row["rev_path"]))
            try:
                for feature in unique_features.itertuples(index=False):
                    score = score_feature(pd.Series(feature._asdict()), bws, sample_cfg)
                    ann = annotation_lookup.get((score["feature_key"], score["sample_id"]))
                    if ann is not None:
                        score.update(
                            {
                                "source_tools": ann.source_tools,
                                "classes": ann.classes,
                                "tool_rows": int(ann.tool_rows),
                            }
                        )
                    records.append(score)
            finally:
                bws.close()

            scored = pd.DataFrame(records)
            with gzip.open(out_path, "wt") as handle:
                scored.to_csv(handle, sep="\t", index=False)

            summary = {
                "sample_id": sample_id,
                "unique_features": int(len(unique_features)),
                "input_rows": int(len(sample_orfs)),
                "scored_features": int((scored["extract_status"] == "ok").sum()) if not scored.empty else 0,
                "coverage_pass": int(scored["coverage_pass"].sum()) if "coverage_pass" in scored else 0,
                "psite_offset": int(sample_cfg.psite_offset),
                "output": str(out_path),
            }
            done.write_text(json.dumps(summary, indent=2) + "\n")
            summaries.append(summary)
            total_unique_scored += int(len(unique_features))
    finally:
        con.close()

    result_paths = sorted(batch_dir.glob("*.scores.tsv.gz"))
    scored_all = pd.concat((pd.read_csv(path, sep="\t") for path in result_paths), ignore_index=True) if result_paths else pd.DataFrame()
    summary = {
        "config": asdict(cfg),
        "db": str(db),
        "bigwig_root": str(bigwig_root),
        "batch_dir": str(batch_dir),
        "samples": summaries,
        "result_rows": int(len(scored_all)),
    }
    (out_dir / "signal_score_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not scored_all.empty:
        scored_all.to_csv(out_dir / "translon_signal_scores.tsv.gz", sep="\t", index=False)
    return scored_all, summary


def expand_scores_by_tool(scored: pd.DataFrame) -> pd.DataFrame:
    """Expand unique feature/sample rows to a tool-level table for plots."""

    if scored.empty or "source_tools" not in scored:
        return scored.copy()
    rows = []
    for row in scored.to_dict("records"):
        tools = str(row.get("source_tools", "")).split("|")
        classes = str(row.get("classes", "")).split("|")
        cls = classes[0] if classes else ""
        for tool in [tool for tool in tools if tool and tool != "nan"]:
            out = dict(row)
            out["source_tool"] = tool
            out["cls"] = cls
            rows.append(out)
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--bigwig-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("figures_bigwig"))
    parser.add_argument("--flank-nt", type=int, default=30)
    parser.add_argument("--body-edge-nt", type=int, default=15)
    parser.add_argument("--min-body-nt", type=int, default=30)
    parser.add_argument("--min-total-signal", type=float, default=10.0)
    parser.add_argument("--max-features-per-sample", type=int, default=2000)
    parser.add_argument("--max-total-features", type=int)
    parser.add_argument("--psite-offsets", type=Path, default=Path("psite_offsets.json"))
    parser.add_argument("--full", action="store_true", help="score all features per sample")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    cfg = ScoreConfig(
        flank_nt=args.flank_nt,
        body_edge_nt=args.body_edge_nt,
        min_body_nt=args.min_body_nt,
        min_total_signal=args.min_total_signal,
        max_features_per_sample=None if args.full else args.max_features_per_sample,
        max_total_features=args.max_total_features,
        overwrite=args.overwrite,
    )
    psite_offsets = load_psite_offsets(args.psite_offsets)
    scored, summary = score_database(args.db, args.bigwig_root, args.out_dir, cfg, psite_offsets=psite_offsets)
    print(json.dumps({k: v for k, v in summary.items() if k != "samples"}, indent=2))
    print(f"Scored rows: {len(scored):,}")


if __name__ == "__main__":
    main()
