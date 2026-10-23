#!/usr/bin/env python3
import argparse
import csv
import shutil
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from read_input_classification import validate_approved


def main():
    p = argparse.ArgumentParser()
    p.add_argument("manifest")
    p.add_argument("validated_output")
    p.add_argument("report")
    a = p.parse_args()
    with Path(a.manifest).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"Approved manifest is empty: {a.manifest}")
    validate_approved(rows)
    source = Path(a.manifest).resolve()
    target = Path(a.validated_output).resolve()
    if source != target:
        shutil.copyfile(source, target)
    Path(a.report).write_text(f"status\tdetail\nok\tvalidated {len(rows)} approved rows\n")


if __name__ == "__main__":
    main()
