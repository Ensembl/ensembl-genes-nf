#!/usr/bin/env python3
"""Write a small, format-aware audit for a native model product."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def count_models(path: Path, native_format: str) -> int:
    if native_format == "bed12":
        return sum(1 for line in path.read_text().splitlines() if line.strip() and not line.startswith("track"))
    transcripts: set[str] = set()
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2] != "exon":
            continue
        for item in fields[8].split(";"):
            key, _, value = item.strip().partition(" ")
            if key == "transcript_id":
                transcripts.add(value.strip('"'))
                break
    return len(transcripts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("native_model", type=Path)
    parser.add_argument("native_format", choices=("bed12", "gtf"))
    parser.add_argument("backend")
    parser.add_argument("scope")
    parser.add_argument("manifest", type=Path)
    parser.add_argument("stats", type=Path)
    parser.add_argument("checksum", type=Path)
    args = parser.parse_args()

    if not args.native_model.is_file() or args.native_model.stat().st_size == 0:
        raise SystemExit(f"native model is missing or empty: {args.native_model}")
    models = count_models(args.native_model, args.native_format)
    if models == 0:
        raise SystemExit(f"native model contains no models: {args.native_model}")
    digest = hashlib.sha256(args.native_model.read_bytes()).hexdigest()
    args.manifest.write_text(
        "backend\tscope\tnative_format\tmodel_path\tmodel_sha256\tstatus\n"
        f"{args.backend}\t{args.scope}\t{args.native_format}\t{args.native_model.name}\t{digest}\tCOMPLETE\n"
    )
    args.stats.write_text(
        "backend\tscope\tnative_format\tmodel_count\tstatus\n"
        f"{args.backend}\t{args.scope}\t{args.native_format}\t{models}\tCOMPLETE\n"
    )
    args.checksum.write_text(f"{digest}  {args.native_model.name}\n")


if __name__ == "__main__":
    main()
