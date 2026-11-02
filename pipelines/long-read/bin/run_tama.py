#!/usr/bin/env python3
"""Run TAMA while making empty or unrepresentable shards explicit."""
from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def has_usable_alignment(path: Path) -> bool:
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            try:
                flag = int(fields[1])
            except ValueError:
                continue
            if flag & 0x4 or fields[2] == "*" or fields[5] == "*":
                continue
            return True
    return False


def is_known_empty_locus_failure(result: subprocess.CompletedProcess[str]) -> bool:
    output = f"{result.stdout}\n{result.stderr}"
    return result.returncode != 0 and "IndexError: list index out of range" in output


def write_empty(path: Path, reason: str) -> None:
    path.write_text("")
    print(f"TAMA produced no models: {reason}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--args", default="")
    parser.add_argument("--tama-executable", default="tama_collapse.py")
    args = parser.parse_args()

    if not has_usable_alignment(args.sam):
        write_empty(args.output, "the shard SAM contains no usable mapped alignments")
        return

    command = [
        args.tama_executable,
        "-s", str(args.sam),
        "-f", str(args.reference),
        "-p", args.prefix,
        "-m", "10",
        *shlex.split(args.args),
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    if result.returncode != 0:
        if is_known_empty_locus_failure(result):
            write_empty(args.output, "TAMA encountered an empty transcript locus")
            return
        raise SystemExit(result.returncode)

    generated = Path(f"{args.prefix}.bed")
    if generated == args.output:
        if not generated.exists():
            write_empty(args.output, "TAMA completed without producing a BED12 file")
    elif generated.exists():
        shutil.copyfile(generated, args.output)
    else:
        write_empty(args.output, "TAMA completed without producing a BED12 file")


if __name__ == "__main__":
    main()
