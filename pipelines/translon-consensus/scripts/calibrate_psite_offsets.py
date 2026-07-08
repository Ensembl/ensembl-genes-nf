#!/usr/bin/env python3
"""Calibrate per-sample mod-3 P-site offsets from MANE_Select CDS signal."""

from __future__ import annotations

import argparse
import json
import logging
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from transcode_bigwig_signal_scores import BigWigPair, discover_bigwigs


ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attrs(attr_text: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attr_text))


def load_mane_cds(gtf: Path) -> list[dict[str, object]]:
    transcripts: dict[str, dict[str, object]] = defaultdict(
        lambda: {"chrom": None, "strand": None, "cds": [], "has_start": False}
    )
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#") or "MANE_Select" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"CDS", "start_codon"}:
                continue
            attrs = parse_attrs(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            record = transcripts[transcript_id]
            chrom = fields[0]
            strand = fields[6]
            if record["chrom"] is None:
                record["chrom"] = chrom
                record["strand"] = strand
            elif record["chrom"] != chrom or record["strand"] != strand:
                logging.warning("Skipping %s with inconsistent chrom/strand", transcript_id)
                transcripts.pop(transcript_id, None)
                continue

            start = int(fields[3]) - 1
            end = int(fields[4])
            if fields[2] == "CDS":
                record["cds"].append((start, end))
            else:
                record["has_start"] = True

    rows = []
    for transcript_id, record in transcripts.items():
        cds = sorted(record["cds"])
        if cds and record["has_start"] and record["chrom"] and record["strand"] in {"+", "-"}:
            rows.append(
                {
                    "transcript_id": transcript_id,
                    "chrom": str(record["chrom"]),
                    "strand": str(record["strand"]),
                    "cds": cds,
                }
            )
    return rows


def load_bigwig_manifest(args: argparse.Namespace) -> pd.DataFrame:
    if args.bigwig_manifest:
        manifest = pd.read_csv(args.bigwig_manifest, sep=None, engine="python")
        required = {"sample_id", "fwd_path", "rev_path"}
        missing = required - set(manifest.columns)
        if missing:
            raise ValueError(f"{args.bigwig_manifest} is missing columns: {', '.join(sorted(missing))}")
        manifest = manifest.copy()
        manifest["has_pair"] = manifest["fwd_path"].notna() & manifest["rev_path"].notna()
    else:
        manifest = discover_bigwigs(args.bigwig_root)

    manifest = manifest[manifest["has_pair"]].copy()
    if args.samples:
        manifest = manifest[manifest["sample_id"].astype(str).isin(set(args.samples))]
    if manifest.empty:
        raise FileNotFoundError("No paired BigWigs found for calibration")
    return manifest.sort_values("sample_id")


def frame_fraction(frames: np.ndarray, frame: int | None = None) -> float:
    total = float(frames.sum())
    if total <= 0:
        return float("nan")
    idx = int(np.argmax(frames)) if frame is None else int(frame)
    return float(frames[idx] / total)


def dominant_frame(frames: np.ndarray) -> int | None:
    if float(frames.sum()) <= 0:
        return None
    return int(np.argmax(frames))


def calibrate_sample(sample_id: str, fwd_path: str, rev_path: str, transcripts: list[dict[str, object]]) -> dict[str, object]:
    frames = {"+": np.zeros(3, dtype=float), "-": np.zeros(3, dtype=float)}
    n_tx = 0
    bws = BigWigPair(fwd_path, rev_path)
    try:
        for tx in transcripts:
            strand = str(tx["strand"])
            arr = bws.transcript_values(str(tx["chrom"]), tx["cds"], strand)
            if arr is None or len(arr) < 60:
                continue
            body = arr[15 : len(arr) - 15]
            if len(body) == 0 or float(body.sum()) <= 0:
                continue
            positions = np.arange(len(body)) % 3
            for frame in range(3):
                frames[strand][frame] += float(body[positions == frame].sum())
            n_tx += 1
    finally:
        bws.close()

    plus_frame = dominant_frame(frames["+"])
    minus_frame = dominant_frame(frames["-"])
    if plus_frame is not None and minus_frame is not None and plus_frame != minus_frame:
        msg = (
            f"{sample_id}: plus and minus dominant frames differ "
            f"({plus_frame} vs {minus_frame}); single-scalar offset is invalid"
        )
        logging.warning(msg)
        raise RuntimeError(msg)

    combined = frames["+"] + frames["-"]
    dom = dominant_frame(combined)
    if dom is None:
        raise RuntimeError(f"{sample_id}: no usable CDS signal for calibration")

    dom_frac = frame_fraction(combined, dom)
    offset = (3 - dom) % 3
    if dom_frac < 0.5:
        logging.warning(
            "%s: dominant frame fraction %.3f is below 0.5; this bigWig is poorly registered",
            sample_id,
            dom_frac,
        )

    return {
        "offset": int(offset),
        "dom_frac": float(dom_frac),
        "dom_frac_plus": frame_fraction(frames["+"], plus_frame),
        "dom_frac_minus": frame_fraction(frames["-"], minus_frame),
        "dom_frame_plus": -1 if plus_frame is None else int(plus_frame),
        "dom_frame_minus": -1 if minus_frame is None else int(minus_frame),
        "n_tx": int(n_tx),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--bigwig-root", type=Path, default=Path("."))
    parser.add_argument("--bigwig-manifest", type=Path)
    parser.add_argument("--samples", nargs="*")
    parser.add_argument("--out", type=Path, default=Path("psite_offsets.json"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    transcripts = load_mane_cds(args.gtf)
    if not transcripts:
        raise RuntimeError(f"No MANE_Select CDS transcripts found in {args.gtf}")

    manifest = load_bigwig_manifest(args)
    results = {}
    print("\t".join(["sample", "dom_frac", "plus_frac", "minus_frac", "plus_frame", "minus_frame", "offset", "n_tx"]))
    for row in manifest.itertuples(index=False):
        sample_id = str(row.sample_id)
        result = calibrate_sample(sample_id, str(row.fwd_path), str(row.rev_path), transcripts)
        results[sample_id] = result
        print(
            "\t".join(
                [
                    sample_id,
                    f"{result['dom_frac']:.3f}",
                    f"{result['dom_frac_plus']:.3f}",
                    f"{result['dom_frac_minus']:.3f}",
                    str(result["dom_frame_plus"]),
                    str(result["dom_frame_minus"]),
                    str(result["offset"]),
                    str(result["n_tx"]),
                ]
            )
        )

    args.out.write_text(json.dumps(results, indent=2) + "\n")
    logging.info("Wrote %s", args.out)


if __name__ == "__main__":
    main()
