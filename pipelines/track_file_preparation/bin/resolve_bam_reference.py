#!/usr/bin/env python3
"""Validate BAM/reference sequence dictionaries and resolve chromosome names."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import List, Sequence, Tuple


def read_bam_dictionary(bam: Path) -> List[Tuple[str, int]]:
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam)], check=True, text=True,
        stdout=subprocess.PIPE,
    )
    rows: List[Tuple[str, int]] = []
    seen = set()
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if not fields or fields[0] != "@SQ":
            continue
        values = {}
        for field in fields[1:]:
            if ":" in field:
                key, value = field.split(":", 1)
                values[key] = value
        name = values.get("SN", "")
        length = values.get("LN", "")
        if not name or not length.isdigit() or int(length) < 1:
            raise ValueError(f"Invalid BAM sequence dictionary row: {line}")
        if name in seen:
            raise ValueError(f"Duplicate BAM sequence name: {name}")
        seen.add(name)
        rows.append((name, int(length)))
    if not rows:
        raise ValueError("BAM sequence dictionary is empty")
    return rows


def read_chrom_sizes(path: Path) -> List[Tuple[str, int]]:
    rows: List[Tuple[str, int]] = []
    seen = set()
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if (len(fields) != 2 or not fields[0] or not fields[1].isdigit()
                    or int(fields[1]) < 1):
                raise ValueError(f"Invalid chrom.sizes row: {line.rstrip()}")
            name, length = fields[0], int(fields[1])
            if name in seen:
                raise ValueError(f"Duplicate chromosome in chrom.sizes: {name}")
            seen.add(name)
            rows.append((name, length))
    if not rows:
        raise ValueError("chrom.sizes is empty")
    return rows


def write_tsv(path: Path, rows: Sequence[Sequence[object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--chrom-sizes", required=True, type=Path)
    parser.add_argument("--gca", required=True)
    parser.add_argument("--assembly-release", required=True)
    parser.add_argument("--entity-id", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()
    prefix = Path(args.output_prefix)
    try:
        bam_rows = read_bam_dictionary(args.bam)
        ref_rows = read_chrom_sizes(args.chrom_sizes)
        bam_by_name = {name: length for name, length in bam_rows}
        ref_by_name = {name: length for name, length in ref_rows}

        if bam_by_name == ref_by_name:
            resolved = sorted(ref_rows)
            resolution = "input_names_match"
        elif set(bam_by_name) == set(ref_by_name):
            raise ValueError("BAM and chrom.sizes names match but sequence lengths differ.")
        elif sorted(bam_by_name.values()) == sorted(ref_by_name.values()):
            resolved = sorted(bam_rows)
            resolution = "bam_names_used_lengths_match"
            print("BAM sequence names differ from chrom.sizes; using BAM names because all sequence lengths match.",
                  file=sys.stderr)
        else:
            raise ValueError("BAM sequence lengths do not match chrom.sizes; refusing to create tracks against a different reference.")

        write_tsv(Path(f"{prefix}.chrom.sizes"), resolved)
        write_tsv(Path(f"{prefix}.reference_resolution.tsv"), [
            ["gca_accession", "assembly_release", "entity_id", "resolution", "chromosome_count"],
            [args.gca, args.assembly_release, args.entity_id, resolution, len(resolved)],
        ])
        version = subprocess.run(["samtools", "--version"], check=True, text=True,
                                 stdout=subprocess.PIPE).stdout.splitlines()[0].split()
        write_tsv(Path(f"{prefix}.reference_resolution.versions.yml"), [
            [f'"RESOLVE_BAM_REFERENCE":'],
            [f"    samtools: {version[1] if len(version) > 1 else 'unknown'}"],
            [f"    reference_resolution: {resolution}"],
        ])
    except (OSError, subprocess.CalledProcessError, ValueError, IndexError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
