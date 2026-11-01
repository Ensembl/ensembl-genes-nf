#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--file-list", required=True, type=Path)
    parser.add_argument("--checksum", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--backend")
    parser.add_argument("--scope")
    parser.add_argument("--status")
    args = parser.parse_args()
    paths = sorted({str(Path(item).resolve()) for item in args.inputs})
    if not paths:
        raise SystemExit("merge requires at least one input")
    args.file_list.write_text("\n".join(paths) + "\n")
    digest = hashlib.sha256(args.file_list.read_bytes()).hexdigest()
    args.checksum.write_text(f"{digest}  {args.file_list.name}\n")
    if args.manifest:
        if not all((args.backend, args.scope, args.status)):
            raise SystemExit("--manifest requires --backend, --scope, and --status")
        args.manifest.write_text(
            json.dumps(
                {
                    "backend": args.backend,
                    "scope": args.scope,
                    "status": args.status,
                    "inputs": paths,
                    "file_list": args.file_list.name,
                    "file_list_sha256": digest,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
