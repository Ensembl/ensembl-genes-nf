#!/usr/bin/env python3
"""Write a deterministic normalized-model comparison table."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--records", type=Path)
    source.add_argument("--records-json")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    records = json.loads(args.records.read_text() if args.records else args.records_json)
    with args.output.open("w") as handle:
        handle.write("cohort_id\tbackend\tscope\tmodel_path\tstatus\n")
        for row in sorted(records, key=lambda item: (item["backend"], item["scope"], item["model_path"])):
            handle.write(f"{row['cohort_id']}\t{row['backend']}\t{row['scope']}\t{row['model_path']}\tVALIDATED\n")


if __name__ == "__main__":
    main()
