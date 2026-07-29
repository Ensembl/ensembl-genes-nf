#!/usr/bin/env python3
"""Reheader a BAM so its reference names match an INSDC FASTA."""

import argparse
import subprocess
import sys
from pathlib import Path


def reference_lengths(fai: Path) -> dict[str, int]:
    lengths = {}
    with fai.open() as handle:
        for line in handle:
            name, length, *_ = line.rstrip("\n").split("\t")
            lengths[name] = int(length)
    return lengths


def report_rows(report: Path) -> list[dict[str, str]]:
    """Read an NCBI assembly_report.txt using its commented column header."""
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
            fields = line.split("\t")
            if columns is None:
                columns = [
                    "Sequence-Name", "Sequence-Role", "Assigned-Molecule",
                    "Assigned-Molecule-Location/Type", "GenBank-Accn",
                    "Relationship", "RefSeq-Accn", "Assembly-Unit",
                    "Sequence-Length", "UCSC-style-name",
                ]
            rows.append(dict(zip(columns, fields)))
    if not rows:
        raise ValueError(f"No assembly rows found in {report}")
    return rows


def bam_sq_names(bam: Path) -> tuple[list[str], list[str]]:
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam)],
        check=True, capture_output=True, text=True,
    )
    sq_lines = [line for line in result.stdout.splitlines() if line.startswith("@SQ\t")]
    names = []
    for line in sq_lines:
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        if "SN" not in fields or "LN" not in fields:
            raise ValueError(f"Malformed @SQ line: {line}")
        names.append((fields["SN"], fields["LN"]))
    if not names:
        raise ValueError(f"BAM has no @SQ records: {bam}")
    return result.stdout.splitlines(), names


def build_mapping(bam_names, target_lengths, report: Path) -> dict[str, str]:
    target_names = set(target_lengths)
    mapping = {}
    used_targets = set()
    rows = report_rows(report)
    for source, bam_length in bam_names:
        if source in target_names:
            target = source
        else:
            target = None
            for row in rows:
                aliases = {value for value in row.values() if value and value != "na"}
                if source in aliases:
                    candidates = list(dict.fromkeys(value for value in row.values() if value in target_names))
                    if len(candidates) == 1:
                        target = candidates[0]
                    elif len(candidates) > 1:
                        raise ValueError(
                            f"Ambiguous INSDC target for BAM reference {source}: {candidates}"
                        )
                    break
        if target is None:
            raise ValueError(
                f"No INSDC FASTA name for BAM reference {source}; check the assembly/report"
            )
        if int(bam_length) != target_lengths[target]:
            raise ValueError(
                f"Length mismatch for {source} -> {target}: "
                f"BAM={bam_length}, INSDC FASTA={target_lengths[target]}"
            )
        if target in used_targets and mapping.get(source) != target:
            raise ValueError(f"Multiple BAM references map to INSDC name {target}")
        mapping[source] = target
        used_targets.add(target)
    return mapping


def rewrite_header(header_lines: list[str], mapping: dict[str, str]) -> str:
    rewritten = []
    for line in header_lines:
        if line.startswith("@SQ\t"):
            fields = line.split("\t")
            for index, field in enumerate(fields):
                if field.startswith("SN:"):
                    fields[index] = "SN:" + mapping[field[3:]]
                    break
            line = "\t".join(fields)
        rewritten.append(line)
    return "\n".join(rewritten) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--reference-fai", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    target_lengths = reference_lengths(args.reference_fai)
    header_lines, bam_names = bam_sq_names(args.bam)
    mapping = build_mapping(bam_names, target_lengths, args.assembly_report)
    new_header = rewrite_header(header_lines, mapping)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as output:
        subprocess.run(
            ["samtools", "reheader", "-", str(args.bam)],
            input=new_header.encode(), stdout=output, check=True,
        )
    print("Reference mapping:")
    for source, target in mapping.items():
        print(f"  {source} -> {target}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"reheader_bam.py: ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
