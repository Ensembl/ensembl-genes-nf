#!/usr/bin/env python3
"""Copy or download a declared input URI and write its checksum metadata."""
from __future__ import annotations

import argparse
import hashlib
import shutil
import urllib.parse
import urllib.request
from pathlib import Path


def fetch(uri: str, output: Path) -> None:
    parsed = urllib.parse.urlparse(uri)
    if parsed.scheme in {"http", "https", "ftp"}:
        with urllib.request.urlopen(uri) as source, output.open("wb") as target:
            shutil.copyfileobj(source, target)
    else:
        source = Path(urllib.parse.unquote(parsed.path if parsed.scheme == "file" else uri))
        shutil.copyfile(source, output)
    if not output.is_file() or output.stat().st_size == 0:
        raise SystemExit(f"acquired input is missing or empty: {uri}")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--uri", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--index-uri")
    parser.add_argument("--index-output", type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--metadata", required=True, type=Path)
    args = parser.parse_args()
    fetch(args.uri, args.output)
    index_name = ""
    index_digest = ""
    if args.index_uri:
        if args.index_output is None:
            raise SystemExit("--index-output is required with --index-uri")
        fetch(args.index_uri, args.index_output)
        index_name = args.index_output.name
        index_digest = digest(args.index_output)
    args.metadata.write_text(
        "id\tinput\tinput_sha256\tinput_index\tinput_index_sha256\n"
        f"{args.run_id}\t{args.output.name}\t{digest(args.output)}\t{index_name}\t{index_digest}\n"
    )


if __name__ == "__main__":
    main()
