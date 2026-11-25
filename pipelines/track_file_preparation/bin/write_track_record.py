#!/usr/bin/env python3
"""Write the per-track result, provenance, exception, and version records."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path
from typing import Optional, Sequence


RESULT_FIELDS = ["gca_accession", "assembly_release", "entity_id", "entity_type",
                 "track_type", "track_path", "source_path", "status", "sha256",
                 "tool_versions", "normalization_parameters"]


def write_tsv(path: Path, header: Sequence[str], row: Sequence[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerow(row)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--track", required=True, type=Path)
    parser.add_argument("--checksum", required=True, type=Path)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--track-result", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--exceptions", required=True, type=Path)
    parser.add_argument("--versions", required=True, type=Path)
    parser.add_argument("--gca", required=True)
    parser.add_argument("--assembly-release", required=True)
    parser.add_argument("--entity-id", required=True)
    parser.add_argument("--entity-type", required=True)
    parser.add_argument("--track-type", required=True)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--source-path", required=True)
    parser.add_argument("--normalization", default="not_applicable")
    parser.add_argument("--stub", action="store_true")
    args = parser.parse_args(argv)

    artifact = Path(args.artifact)
    artifact.parent.mkdir(parents=True, exist_ok=True)
    args.track_result.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.exceptions.parent.mkdir(parents=True, exist_ok=True)
    args.versions.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.track, artifact)
    checksum = "stub" if args.stub else args.checksum.read_text().split()[0]
    write_tsv(args.track_result, RESULT_FIELDS, [
        args.gca, args.assembly_release, args.entity_id, args.entity_type,
        args.track_type, args.destination, args.source_path, "complete", checksum,
        "tracked", args.normalization,
    ])
    shutil.copyfile(args.track_result, args.provenance)
    args.exceptions.write_text("")
    version = "stub" if args.stub else "pipeline"
    args.versions.write_text(f'"WRITE_TRACK_RECORD":\n    record_writer: {version}\n')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
