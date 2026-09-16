#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("manifest")
    p.add_argument("validated_output")
    p.add_argument("report")
    a = p.parse_args()
    lines = [line for line in Path(a.manifest).read_text().splitlines() if line.strip()]
    if len(lines) <= 1:
        raise SystemExit(f"Approved manifest is empty: {a.manifest}")
    shutil.copyfile(a.manifest, a.validated_output)
    Path(a.report).write_text("status\tdetail\nok\tapproved manifest validated\n")


if __name__ == "__main__":
    main()
