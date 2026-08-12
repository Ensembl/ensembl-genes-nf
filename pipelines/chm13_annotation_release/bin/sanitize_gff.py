#!/usr/bin/env python3
"""Normalize projected and HAVANA GFF3 for the Ensembl core loader."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from urllib.parse import unquote


def parse_attributes(raw: str) -> list[tuple[str, str | None]]:
    result = []
    for item in raw.rstrip("\n").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            result.append((key, value))
        else:
            result.append((item, None))
    return result


def sanitize(input_path: Path, output_path: Path, kind: str) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with input_path.open(encoding="utf-8") as source, output_path.open("w", encoding="utf-8") as target:
        for raw in source:
            if raw.startswith("#"):
                target.write(raw)
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                counts["skipped_malformed"] += 1
                continue
            feature = fields[2]
            attributes = parse_attributes(fields[8])
            values = {key: unquote(value or "") for key, value in attributes}

            if kind == "manual" and feature == "gene_segment":
                fields[2] = "transcript"
                feature = "transcript"
                counts["gene_segment_to_transcript"] += 1

            if kind == "projected" and feature in {"gene", "transcript"}:
                biotype = values.get("biotype")
                source_key = "gene_type" if feature == "gene" else "transcript_type"
                if not biotype and values.get(source_key):
                    attributes.append(("biotype", values[source_key]))
                    counts[f"{source_key}_to_biotype"] += 1

            fields[8] = ";".join(
                key if value is None else f"{key}={value}" for key, value in attributes
            )
            target.write("\t".join(fields) + "\n")
            counts["records"] += 1
            counts[f"feature_{feature}"] += 1

    return dict(counts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kind", choices=("manual", "projected"), required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = sanitize(args.input, args.output, args.kind)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
