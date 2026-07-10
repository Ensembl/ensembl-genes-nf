#!/usr/bin/env python3
"""Prepare annotation and offset inputs required by ChOROS.

The pipeline's RiboMetric annotation is already transcript-relative.  ChOROS
uses the same information with different column names and requires an A-site
offset for every (read length, 5-prime frame) pair.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def nearest_frame_offset(base_offset: int, frame: int, read_length: int) -> int:
    """Return the closest valid offset placing the A-site in CDS frame zero."""
    candidates = [
        offset
        for offset in range(read_length - 2)
        if (frame + offset) % 3 == 0
    ]
    if not candidates:
        raise ValueError(f"no valid A-site offset for read length {read_length}")
    return min(candidates, key=lambda offset: (abs(offset - base_offset), offset))


def read_best_offsets(path: Path) -> dict[int, int]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"length", "offset"}.issubset(reader.fieldnames or []):
            raise ValueError("offset file must contain length and offset columns")
        offsets = {int(row["length"]): int(row["offset"]) for row in reader}
    if not offsets:
        raise ValueError("offset file contains no usable rows")
    return offsets


def write_choros_offsets(offsets: dict[int, int], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["length", "frame_0", "frame_1", "frame_2"])
        for read_length, base_offset in sorted(offsets.items()):
            if not 0 <= base_offset < read_length - 2:
                raise ValueError(
                    f"offset {base_offset} is invalid for read length {read_length}"
                )
            writer.writerow(
                [read_length]
                + [
                    nearest_frame_offset(base_offset, frame, read_length)
                    for frame in range(3)
                ]
            )


def write_choros_lengths(annotation: Path, path: Path) -> None:
    with annotation.open(newline="") as source, path.open("w", newline="") as dest:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"transcript_id", "cds_start", "cds_end", "transcript_length"}
        if not required.issubset(reader.fieldnames or []):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"RiboMetric annotation is missing: {', '.join(missing)}")
        writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript", "utr5_length", "cds_length", "utr3_length"])
        seen = set()
        for row in reader:
            transcript = row["transcript_id"]
            if transcript in seen:
                continue
            seen.add(transcript)
            cds_start = int(row["cds_start"])
            cds_end = int(row["cds_end"])
            transcript_length = int(row["transcript_length"])
            if not 0 <= cds_start < cds_end <= transcript_length:
                raise ValueError(
                    f"invalid transcript-relative CDS coordinates for {transcript}"
                )
            writer.writerow(
                [transcript, cds_start, cds_end - cds_start, transcript_length - cds_end]
            )
        if not seen:
            raise ValueError("RiboMetric annotation contains no transcripts")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--offsets", required=True, type=Path)
    parser.add_argument("--lengths-output", required=True, type=Path)
    parser.add_argument("--offsets-output", required=True, type=Path)
    args = parser.parse_args()

    write_choros_lengths(args.annotation, args.lengths_output)
    write_choros_offsets(read_best_offsets(args.offsets), args.offsets_output)


if __name__ == "__main__":
    main()
