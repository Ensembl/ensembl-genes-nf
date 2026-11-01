#!/usr/bin/env python3
"""Build a FASTA sequence dictionary and content identity record."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def dictionary(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    name: str | None = None
    length = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    rows.append((name, length))
                name = line[1:].split(None, 1)[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        rows.append((name, length))
    if not rows:
        raise SystemExit(f"reference FASTA contains no sequences: {path}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--dictionary", required=True, type=Path)
    parser.add_argument("--identity", required=True, type=Path)
    args = parser.parse_args()
    rows = dictionary(args.reference)
    args.dictionary.write_text("".join(f"{name}\t{length}\n" for name, length in rows))
    args.identity.write_text(json.dumps({
        "path": str(args.reference),
        "sha256": digest(args.reference),
        "dictionary_sha256": digest(args.dictionary),
    }) + "\n")


if __name__ == "__main__":
    main()
