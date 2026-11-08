#!/usr/bin/env python3
"""Prepare a reproducible SLURM/Apptainer layout for translon characterisation.

This is an environment and preflight helper, not a job submitter. It creates
directories plus an explicit Nextflow config/report below the requested root;
it never downloads or overwrites reference data.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


MIN_NEXTFLOW = (26, 4, 6)
REQUIRED_ASSETS = (
    "pipelines/translon-characterisation/assets/production-reference-contract.json",
    "pipelines/translon-characterisation/assets/hg38.multiz100way.primary.maf.jsonl",
    "pipelines/translon-characterisation/assets/hg38.multiz100way.58mammals.tsv",
)
DIRECTORIES = (
    "apptainer-cache",
    "tmp",
    "work",
    "logs",
    "results",
    "inputs",
    "references",
    "models",
)


def command_version(command: str, args: tuple[str, ...] = ("--version",)) -> str | None:
    executable = shutil.which(command)
    if not executable:
        return None
    result = subprocess.run(
        [executable, *args], capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        return None
    return (result.stdout + result.stderr).strip() if (result.stdout + result.stderr).strip() else ""


def parse_nextflow_version(value: str | None) -> tuple[int, ...] | None:
    if not value:
        return None
    match = re.search(r"version\s+(\d+(?:\.\d+)+)", value, re.IGNORECASE)
    if not match:
        return None
    return tuple(int(part) for part in match.group(1).split("."))


def parse_java_major(value: str | None) -> int | None:
    if not value:
        return None
    match = re.search(r'version\s+"(\d+)', value)
    return int(match.group(1)) if match else None


def runtime_checks(skip: bool, executor: str) -> tuple[dict[str, str | None], list[str], list[str]]:
    versions: dict[str, str | None] = {}
    errors: list[str] = []
    warnings: list[str] = []
    if skip:
        return versions, errors, warnings

    versions["nextflow"] = command_version("nextflow", ("-version",))
    nextflow_version = parse_nextflow_version(versions["nextflow"])
    if nextflow_version is None:
        errors.append("nextflow is missing or its version could not be parsed")
    elif nextflow_version < MIN_NEXTFLOW:
        errors.append(
            f"nextflow {'.'.join(map(str, MIN_NEXTFLOW))}+ is required; found {versions['nextflow']}"
        )

    versions["java"] = command_version("java", ("-version",))
    java_major = parse_java_major(versions["java"])
    if java_major is None:
        errors.append("Java is missing or its version could not be parsed")
    elif java_major not in {17, 21}:
        errors.append(f"Java 17 or 21 is required; found Java {java_major}")

    apptainer = shutil.which("apptainer")
    singularity = shutil.which("singularity")
    if apptainer:
        versions["apptainer"] = command_version("apptainer")
        if not versions["apptainer"]:
            errors.append("Apptainer is present but could not execute; start/configure it before submission")
    elif singularity:
        versions["singularity"] = command_version("singularity")
        if not versions["singularity"]:
            errors.append("Singularity is present but could not execute; start/configure it before submission")
        warnings.append("Apptainer was not found; Singularity will be used instead")
    else:
        errors.append("Apptainer or Singularity is required for the slurm_apptainer profile")

    versions["sbatch"] = shutil.which("sbatch")
    if executor == "slurm" and not versions["sbatch"]:
        warnings.append("sbatch was not found; this host cannot submit SLURM jobs")
    return versions, errors, warnings


def validate_repo(repo: Path) -> list[str]:
    missing = []
    for relative in REQUIRED_ASSETS:
        if not (repo / relative).is_file():
            missing.append(str(repo / relative))
    return missing


def validate_manifest(path: Path | None) -> tuple[list[dict[str, object]], list[str]]:
    if path is None:
        return [], []
    if not path.is_file():
        return [], [f"input manifest does not exist: {path}"]
    records: list[dict[str, object]] = []
    errors: list[str] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"name", "path"}
        if not required_columns <= set(reader.fieldnames or ()):
            return [], [f"input manifest must contain columns: {', '.join(sorted(required_columns))}"]
        for row in reader:
            if not row.get("name") or not row.get("path"):
                continue
            source = Path(row["path"]).expanduser().resolve()
            required = str(row.get("required", "true")).lower() not in {"0", "false", "no"}
            record = {
                "name": row["name"],
                "path": str(source),
                "required": required,
                "exists": source.exists(),
                "size": source.stat().st_size if source.is_file() else None,
            }
            records.append(record)
            if required and not source.is_file():
                errors.append(f"required input is missing: {source}")
    return records, errors


def write_if_allowed(path: Path, content: str, force: bool) -> bool:
    if path.exists() and not force:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return True


def hpc_config_text(root: Path, executor: str, partition: str | None, account: str | None) -> str:
    lines = [
        "// Generated by setup_hpc_environment.py; pass with -c on the Nextflow command line.",
        "process {",
        f"    executor = {executor!r}",
    ]
    if executor == "slurm" and partition:
        lines.append(f"    queue = {partition!r}")
    if executor == "slurm" and account:
        lines.append(f"    clusterOptions = '--account={account}'")
    lines.extend([
        "}",
        "apptainer {",
        "    enabled = true",
        "    autoMounts = true",
        f"    cacheDir = {str(root / 'apptainer-cache')!r}",
        "    runOptions = '--cleanenv'",
        "}",
        "singularity.enabled = false",
    ])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, help="Shared no-backup root for this run")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3], help="ensembl-genes-nf checkout")
    parser.add_argument("--executor", choices=("local", "slurm"), default="slurm",
                        help="Run tasks inside an existing salloc ('local') or submit child jobs ('slurm')")
    parser.add_argument("--partition", help="SLURM partition/queue")
    parser.add_argument("--account", help="SLURM account passed as --account")
    parser.add_argument("--inputs-manifest", type=Path, help="Optional TSV with name, path, and required columns")
    parser.add_argument("--check-only", action="store_true", help="Validate tools and inputs without creating files")
    parser.add_argument("--skip-runtime-check", action="store_true", help="Skip Nextflow/Java/Apptainer checks")
    parser.add_argument("--force", action="store_true", help="Overwrite generated hpc.config and report")
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    repo = args.repo.expanduser().resolve()

    versions, errors, warnings = runtime_checks(args.skip_runtime_check, args.executor)
    errors.extend(f"required repository asset is missing: {path}" for path in validate_repo(repo))
    manifest_records, manifest_errors = validate_manifest(args.inputs_manifest)
    errors.extend(manifest_errors)

    report = {
        "root": str(root),
        "repo": str(repo),
        "profile": "apptainer_local" if args.executor == "local" else "slurm_apptainer",
        "executor": args.executor,
        "config": str(root / "hpc.config"),
        "versions": versions,
        "directories": [str(root / relative) for relative in DIRECTORIES],
        "input_manifest": manifest_records,
        "errors": errors,
        "warnings": warnings,
    }
    if errors:
        print(json.dumps(report, indent=2, sort_keys=True), file=sys.stderr)
        return 2
    if args.check_only:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    for relative in DIRECTORIES:
        (root / relative).mkdir(parents=True, exist_ok=True)
    config_path = root / "hpc.config"
    report_path = root / "setup-report.json"
    write_if_allowed(config_path, hpc_config_text(root, args.executor, args.partition, args.account), args.force)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(f"Prepared {root}")
    print(f"Use this explicit Nextflow config with: -c {config_path}")
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
