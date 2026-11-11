#!/usr/bin/env python3
"""Fetch or reuse one assembly's reference files in a persistent cache."""

import argparse
import gzip
import hashlib
import shutil
import urllib.request
from pathlib import Path


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def expected_digest(value):
    if not value or value == "-":
        return ""
    return value.removeprefix("sha256:").strip().lower()


def fetch(source, target, checksum, decompress=False):
    if not source or source == "-":
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    marker = target.with_suffix(target.suffix + ".source")
    checksum = expected_digest(checksum)
    marker_value = f"{source}\t{checksum}\n"
    if target.exists() and marker.exists() and marker.read_text() == marker_value:
        return

    temporary = target.with_name(target.name + ".download")
    if source.startswith(("http://", "https://", "ftp://")):
        with urllib.request.urlopen(source) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
    else:
        source_path = Path(source).expanduser().resolve()
        if not source_path.is_file():
            raise SystemExit(f"Reference source does not exist: {source}")
        shutil.copyfile(source_path, temporary)

    actual = digest(temporary)
    if checksum and actual != checksum:
        temporary.unlink(missing_ok=True)
        raise SystemExit(f"Checksum mismatch for {source}: expected {checksum}, got {actual}")
    if decompress or source.lower().split("?", 1)[0].endswith(".gz"):
        decompressed = target.with_name(target.name + ".partial")
        with gzip.open(temporary, "rb") as source_handle, decompressed.open("wb") as target_handle:
            shutil.copyfileobj(source_handle, target_handle, length=1024 * 1024)
        temporary.unlink()
        decompressed.replace(target)
    else:
        temporary.replace(target)
    marker.write_text(marker_value)


def derive_chrom_sizes(fasta, target):
    marker = target.with_suffix(target.suffix + ".source")
    marker_value = f"derived-from-fasta\t{digest(fasta)}\n"
    if target.exists() and marker.exists() and marker.read_text() == marker_value:
        return
    sizes = []
    name = None
    length = 0
    with fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sizes.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        sizes.append((name, length))
    target.write_text("".join(f"{name}\t{length}\n" for name, length in sizes))
    marker.write_text(marker_value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--fasta-source", required=True)
    parser.add_argument("--gtf-source", required=True)
    parser.add_argument("--chrom-sizes-source", required=True)
    parser.add_argument("--fasta-checksum", default="")
    parser.add_argument("--gtf-checksum", default="")
    parser.add_argument("--chrom-sizes-checksum", default="")
    args = parser.parse_args()

    cache = args.cache_dir / args.assembly
    fetch(args.fasta_source, cache / "reference.fa", args.fasta_checksum)
    fetch(args.gtf_source, cache / "reference.gtf", args.gtf_checksum)
    if args.chrom_sizes_source and args.chrom_sizes_source != "-":
        fetch(args.chrom_sizes_source, cache / "chrom.sizes", args.chrom_sizes_checksum)
    else:
        derive_chrom_sizes(cache / "reference.fa", cache / "chrom.sizes")
    shutil.copyfile(cache / "reference.fa", "reference.fa")
    shutil.copyfile(cache / "reference.gtf", "reference.gtf")
    shutil.copyfile(cache / "chrom.sizes", "chrom.sizes")


if __name__ == "__main__":
    main()
