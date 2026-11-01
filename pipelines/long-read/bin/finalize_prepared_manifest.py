#!/usr/bin/env python3
"""Replace unresolved/source manifest paths with checksum-bearing prepared artifacts."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--normalized", required=True, type=Path)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--checksums", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    args = parser.parse_args()

    with args.normalized.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    prepared = {}
    with args.prepared.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("id") == "id":
                continue
            prepared[row["id"]] = row
    if not rows:
        raise SystemExit("normalized manifest contains no rows")
    missing = [row["id"] for row in rows if row["id"] not in prepared]
    if missing:
        raise SystemExit("prepared inputs are missing rows: " + ", ".join(missing))

    fieldnames = list(rows[0])
    for name in ("input", "input_sha256", "input_index", "input_index_sha256"):
        if name not in fieldnames:
            fieldnames.append(name)
    output_rows = []
    args.checksums.write_text("id\tpath\tsha256\n")
    with args.checksums.open("a") as checksums:
        for row in rows:
            update = prepared[row["id"]]
            row["input"] = update["input"]
            row["input_sha256"] = update["input_sha256"]
            row["input_index"] = update.get("input_index", "")
            row["input_index_sha256"] = update.get("input_index_sha256", "")
            checksums.write(f"{row['id']}\t{row['input']}\t{row['input_sha256']}\n")
            if row["input_index"]:
                checksums.write(f"{row['id']}[index]\t{row['input_index']}\t{row['input_index_sha256']}\n")
            output_rows.append(row)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)
    args.output.with_suffix(".sha256").write_text(f"{sha256(args.output)}  {args.output.name}\n")
    with args.status.open("w") as handle:
        handle.write("id\tstatus\tmessage\n")
        for row in output_rows:
            handle.write(f"{row['id']}\tVALID\tprepared input checksum recorded\n")


if __name__ == "__main__":
    main()
