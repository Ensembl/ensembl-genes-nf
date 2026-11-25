#!/usr/bin/env python3
"""Validate a STAR SJ.out.tab (or the documented nine-column TSV) file."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


HEADER = ["chrom", "start", "end", "strand", "motif", "annotated",
          "uniqReads", "multiReads", "maxOverhang"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--junctions", required=True, type=Path)
    parser.add_argument("--entity-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        rows = 0
        with args.junctions.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for line_number, fields in enumerate(reader, 1):
                if line_number == 1 and fields == HEADER:
                    continue
                rows += 1
                try:
                    start, end = int(fields[1]), int(fields[2])
                except (IndexError, ValueError):
                    raise ValueError(f"Malformed STAR junction row: {' '.join(fields)}")
                if len(fields) != 9 or start < 1 or end < start:
                    raise ValueError(f"Malformed STAR junction row: {' '.join(fields)}")
        if rows == 0:
            raise ValueError("STAR junction file has no data rows")
        args.output.write_text(
            "entity_id\ttrack_type\tstatus\n"
            f"{args.entity_id}\tsplice_junction\tcomplete\n"
        )
    except (OSError, ValueError, csv.Error) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
