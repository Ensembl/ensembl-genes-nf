#!/usr/bin/env python3
"""Build a replacement BAM header using an INSDC FASTA index and assembly report."""

import argparse
import re
import sys
from pathlib import Path
from typing import Optional


def reference_lengths(fai: Path) -> dict[str, int]:
    lengths = {}
    with fai.open() as handle:
        for line in handle:
            name, length, *_ = line.rstrip("\n").split("\t")
            lengths[name] = int(length)
    return lengths


def report_rows(report: Path) -> list[dict[str, str]]:
    columns = None
    rows = []
    with report.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                candidate = line.lstrip("# ").split("\t")
                if "Sequence-Name" in candidate:
                    columns = candidate
                continue
            if not line.strip():
                continue
            if columns is None:
                columns = [
                    "Sequence-Name", "Sequence-Role", "Assigned-Molecule",
                    "Assigned-Molecule-Location/Type", "GenBank-Accn",
                    "Relationship", "RefSeq-Accn", "Assembly-Unit",
                    "Sequence-Length", "UCSC-style-name",
                ]
            rows.append(dict(zip(columns, line.split("\t"))))
    if not rows:
        raise ValueError(f"No assembly rows found in {report}")
    return rows


def header_lines(header: Path) -> tuple[list[str], list[tuple[str, int]]]:
    lines = header.read_text().splitlines()
    names = []
    for line in lines:
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        if "SN" not in fields or "LN" not in fields:
            raise ValueError(f"Malformed @SQ line: {line}")
        names.append((fields["SN"], int(fields["LN"])))
    if not names:
        raise ValueError(f"Header has no @SQ records: {header}")
    return lines, names


def canonical_name(value: str) -> str:
    """Normalise common Ensembl/INSDC chromosome aliases for lookup only."""
    value = value.strip().lower()
    value = re.sub(r"^chr(?=(?:\d+|x|y|m|mt|un))", "", value)
    if value in {"m", "mt", "mitochondrion", "mitochondrial", "mitochondrial_genome"}:
        return "__mitochondrial__"
    return value


def row_target_names(row: dict[str, str], target_names: set[str]) -> list[str]:
    """Return target FASTA names represented by one assembly-report row."""
    candidates = list(dict.fromkeys(value for value in row.values() if value in target_names))
    target_by_canonical = {}
    for target in target_names:
        target_by_canonical.setdefault(canonical_name(target), []).append(target)
    for value in row.values():
        for target in target_by_canonical.get(canonical_name(value), []):
            if target not in candidates:
                candidates.append(target)
    return candidates


def fasta_name_rows(names_file: Path) -> list[dict[str, str]]:
    rows = []
    with names_file.open() as handle:
        for line in handle:
            target, _, description = line.rstrip("\n").partition("\t")
            if target:
                tokens = re.findall(r"[A-Za-z0-9_.-]+", description)
                rows.append({"Reference-Name": target, "Description": description, **{
                    f"Description-{index}": token for index, token in enumerate(tokens)
                }})
    return rows


def build_mapping(bam_names, target_lengths, report: Path, names_file: Optional[Path] = None) -> dict[str, str]:
    target_names = set(target_lengths)
    mapping = {}
    used_targets = set()
    rows = report_rows(report)
    if names_file:
        rows.extend(fasta_name_rows(names_file))
    for source, bam_length in bam_names:
        if source in target_names:
            target = source
        else:
            target = None
            matching_rows = [
                row for row in rows
                if source in row.values()
                or canonical_name(source) in {canonical_name(value) for value in row.values()}
            ]
            for row in matching_rows:
                candidates = row_target_names(row, target_names)
                if len(candidates) == 1:
                    target = candidates[0]
                elif len(candidates) > 1:
                    raise ValueError(f"Ambiguous INSDC target for {source}: {candidates}")
                if target is not None:
                    break
        if target is None:
            raise ValueError(f"No INSDC FASTA name for {source}; check assembly/report")
        if bam_length != target_lengths[target]:
            raise ValueError(
                f"Length mismatch for {source} -> {target}: "
                f"BAM={bam_length}, INSDC FASTA={target_lengths[target]}"
            )
        if target in used_targets and mapping.get(source) != target:
            raise ValueError(f"Multiple BAM references map to INSDC name {target}")
        mapping[source] = target
        used_targets.add(target)
    return mapping


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--header", type=Path, required=True)
    parser.add_argument("--reference-fai", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--reference-names", type=Path, required=False)
    parser.add_argument("--output-header", type=Path, required=True)
    args = parser.parse_args()

    lines, bam_names = header_lines(args.header)
    mapping = build_mapping(
        bam_names, reference_lengths(args.reference_fai), args.assembly_report, args.reference_names
    )
    rewritten = []
    for line in lines:
        if line.startswith("@SQ\t"):
            fields = line.split("\t")
            for index, field in enumerate(fields):
                if field.startswith("SN:"):
                    fields[index] = "SN:" + mapping[field[3:]]
                    break
            line = "\t".join(fields)
        rewritten.append(line)
    args.output_header.write_text("\n".join(rewritten) + "\n")
    for source, target in mapping.items():
        print(f"{source} -> {target}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"build_reheader_header.py: ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
