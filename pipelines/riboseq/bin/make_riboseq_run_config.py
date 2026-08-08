#!/usr/bin/env python3
"""Create a validated Ribo-seq run-context Nextflow config.

The generated config contains only stable input/output locations and reference
paths. Run-specific processing choices can still be supplied on the command
line, so the file is safe to reuse for retries and fresh runs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path


def find_one(root: Path, names: list[str], patterns: list[str], label: str) -> Path:
    candidates = [root / name for name in names]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    for pattern in patterns:
        matches = sorted(path for path in root.glob(pattern) if path.exists())
        if len(matches) == 1:
            return matches[0].resolve()
        if len(matches) > 1:
            raise ValueError(
                f"Multiple {label} candidates found under {root}: "
                + ", ".join(str(path) for path in matches)
            )

    raise ValueError(f"Could not find {label} under {root}")


def validate_sample_sheet(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix not in {".csv", ".tsv"}:
        raise ValueError(f"Sample sheet must have .csv or .tsv suffix: {path}")

    delimiter = "\t" if suffix == ".tsv" else ","
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        required = {"Run", "study_accession"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Sample sheet {path} is missing required columns: "
                + ", ".join(sorted(missing))
            )
        rows = list(reader)

    if not rows:
        raise ValueError(f"Sample sheet is empty: {path}")
    if any(not row["Run"] or not row["study_accession"] for row in rows):
        raise ValueError("Sample sheet contains an empty Run or study_accession value")
    return delimiter


def nf_string(value: Path | str) -> str:
    return json.dumps(str(value))


def read_reference_config(path: Path) -> dict[str, Path]:
    """Read the path-valued reference parameters from a generated config."""
    required = {
        "star_index",
        "gtf",
        "fasta",
        "chrom_sizes_file",
        "ribometric_annotation",
        "transcriptome_fasta",
    }
    values: dict[str, Path] = {}
    assignment = re.compile(r"^\s*(star_index|gtf|fasta|chrom_sizes_file|ribometric_annotation|transcriptome_fasta)\s*=\s*(['\"])(.*?)\2\s*$")
    for line in path.read_text().splitlines():
        match = assignment.match(line)
        if match:
            value = Path(match.group(3))
            values[match.group(1)] = (path.parent / value).resolve() if not value.is_absolute() else value.resolve()

    missing = required - values.keys()
    if missing:
        raise ValueError(
            f"Reference config {path} is missing required path parameters: "
            + ", ".join(sorted(missing))
        )
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", required=True, type=Path)
    parser.add_argument(
        "--reference-config",
        type=Path,
        help="Existing generated params.config; defaults to params.config under --reference-dir",
    )
    parser.add_argument("--collapsed-read-path", required=True, type=Path)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--output", type=Path, help="Output .config path")
    parser.add_argument("--force", action="store_true", help="Replace an existing config")
    args = parser.parse_args()

    reference_dir = args.reference_dir.resolve()
    collapsed_read_path = args.collapsed_read_path.resolve()
    sample_sheet = args.sample_sheet.resolve()
    outdir = args.outdir.resolve()

    for path, label in [
        (reference_dir, "reference directory"),
        (collapsed_read_path, "collapsed-read directory"),
        (sample_sheet, "sample sheet"),
    ]:
        if not path.exists():
            raise ValueError(f"{label} does not exist: {path}")

    delimiter = validate_sample_sheet(sample_sheet)
    reference_config = args.reference_config
    reference_config_explicit = reference_config is not None
    if reference_config is None:
        for candidate in (reference_dir / "params.config", reference_dir / "riboseq_params.config"):
            if candidate.is_file():
                reference_config = candidate
                break

    if reference_config is not None:
        reference_config = reference_config.resolve()
        if not reference_config.is_file():
            raise ValueError(f"Reference config does not exist: {reference_config}")
        try:
            references = read_reference_config(reference_config)
        except ValueError:
            if reference_config_explicit:
                raise
            references = None
    else:
        references = None

    if references is None:
        references = {
            "star_index": find_one(reference_dir, ["star_index"], ["**/star_index"], "STAR index"),
            "gtf": find_one(reference_dir, ["annotation.gtf", "genes.gtf"], ["**/*.gtf"], "GTF"),
            "fasta": find_one(reference_dir, ["genome.fa", "genome.fasta"], ["**/genome*.fa", "**/genome*.fasta"], "genome FASTA"),
            "chrom_sizes_file": find_one(reference_dir, ["genome.chrom.sizes", "chrom.sizes"], ["**/*.chrom.sizes", "**/chrom.sizes"], "chromosome sizes"),
            "ribometric_annotation": find_one(reference_dir, [], ["**/*ribometric*.tsv", "**/*RiboMetric*.tsv"], "RiboMetric annotation"),
            "transcriptome_fasta": find_one(reference_dir, [], ["**/*transcripts*.fa", "**/*transcriptome*.fa"], "transcriptome FASTA"),
        }

    for key, path in references.items():
        if not path.exists():
            raise ValueError(f"Reference parameter {key} does not exist: {path}")

    output = (args.output or outdir / "pipeline_info" / "run_context.config").resolve()
    if output.exists() and not args.force:
        raise ValueError(f"Output already exists; use --force to replace it: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    config_lines = [
        "// Generated by make_riboseq_run_config.py",
        f"// Sample-sheet delimiter: {'TSV' if delimiter == chr(9) else 'CSV'}",
        "params {",
        f"    sample_sheet = {nf_string(sample_sheet)}",
        f"    sample_sheet_sep = {nf_string(delimiter)}",
        f"    collapsed_read_path = {nf_string(collapsed_read_path)}",
        f"    outdir = {nf_string(outdir)}",
    ]
    config_lines.extend(f"    {key} = {nf_string(value)}" for key, value in references.items())
    config_lines.extend(["}", ""])
    output.write_text("\n".join(config_lines))

    context = {
        "sample_sheet": str(sample_sheet),
        "sample_sheet_delimiter": "tsv" if delimiter == "\t" else "csv",
        "collapsed_read_path": str(collapsed_read_path),
        "outdir": str(outdir),
        "references": {key: str(value) for key, value in references.items()},
    }
    (output.parent / "run_context.json").write_text(json.dumps(context, indent=2) + "\n")
    print(output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
