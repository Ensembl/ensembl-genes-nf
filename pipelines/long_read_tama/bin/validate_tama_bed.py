#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("bed")
    p.add_argument("report")
    p.add_argument("--validated-output")
    p.add_argument("--run")
    p.add_argument("--shard")
    a = p.parse_args()
    source = Path(a.bed)
    if not source.is_file() or source.stat().st_size == 0:
        raise SystemExit(f"TAMA BED is missing or empty: {source}")
    rows = [line for line in source.read_text().splitlines() if line.strip()]
    bad = [i for i, line in enumerate(rows, 1) if len(line.split()) != 12]
    if bad:
        raise SystemExit(f"Invalid TAMA BED12 rows: {len(bad)}")
    if a.validated_output:
        shutil.copyfile(source, a.validated_output)
    if a.run is None:
        Path(a.report).write_text(f"models\t{len(rows)}\n")
    else:
        Path(a.report).write_text(f"run\tshard\tmodels\n{a.run}\t{a.shard}\t{len(rows)}\n")


if __name__ == "__main__":
    main()
