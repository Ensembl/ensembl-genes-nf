#!/usr/bin/env python3
"""Write a deterministic FLAIR input/native-product manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, *, root: Path | None = None) -> dict[str, object]:
    name = str(path.relative_to(root)) if root is not None else path.name
    return {"name": name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--native-products", type=Path)
    parser.add_argument("--run-id", default="")
    parser.add_argument("--shard-id", default="")
    parser.add_argument("--status", required=True)
    parser.add_argument("--status-file", type=Path)
    parser.add_argument("--args", default="")
    parser.add_argument("--input", nargs="*", default=[])
    args = parser.parse_args()

    if args.status_file:
        rows = [line for line in args.status_file.read_text().splitlines() if line.strip()]
        if not rows:
            raise SystemExit(f"status file is empty: {args.status_file}")
        args.status = rows[-1].split("\t")[-1]

    inputs = []
    for value in sorted(args.input):
        path = Path(value)
        if path.is_file():
            inputs.append(file_record(path))

    products = []
    if args.native_products and args.native_products.is_dir():
        for path in sorted(p for p in args.native_products.rglob("*") if p.is_file()):
            products.append(file_record(path, root=args.native_products))

    payload = {
        "backend": "flair",
        "run_id": args.run_id,
        "shard_id": args.shard_id,
        "status": args.status,
        "args": args.args,
        "inputs": inputs,
        "native_products": products,
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
