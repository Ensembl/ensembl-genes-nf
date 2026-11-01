#!/usr/bin/env python3
"""Run a long-read pipeline benchmark and retain reproducibility evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import shlex
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def output_checksums(root: Path) -> list[dict[str, object]]:
    records = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        records.append({"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    return records


def trace_summary(trace: Path) -> dict[str, object]:
    if not trace.is_file():
        return {"trace_present": False, "tasks": 0, "retries": 0}
    with trace.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    retries = 0
    for row in rows:
        try:
            retries += int(row.get("retries", 0) or 0)
        except ValueError:
            pass
    resource_fields = ["task_id", "hash", "name", "status", "realtime", "%cpu", "peak_rss", "peak_vmem", "cpus", "memory", "retries"]
    resource_rows = [{field: row.get(field, "") for field in resource_fields} for row in rows]
    return {"trace_present": True, "tasks": len(rows), "retries": retries, "resource_rows": resource_rows}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--params-file", type=Path)
    parser.add_argument("--nextflow", default="nextflow")
    parser.add_argument("--profile")
    parser.add_argument("--extra", nargs=argparse.REMAINDER, default=[])
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    command = [args.nextflow, "run", str(args.pipeline), "--manifest", str(args.manifest), "--reference_fasta", str(args.reference)]
    if args.params_file:
        command.extend(["-params-file", str(args.params_file)])
    if args.profile:
        command.extend(["-profile", args.profile])
    command.extend(["--outdir", str(args.outdir)])
    command.extend(args.extra)

    stdout_path = args.outdir / "benchmark.stdout.log"
    stderr_path = args.outdir / "benchmark.stderr.log"
    started_at = datetime.now(timezone.utc)
    start = time.monotonic()
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, text=True, check=False)
    finished_at = datetime.now(timezone.utc)

    checksums = output_checksums(args.outdir)
    (args.outdir / "benchmark_checksums.tsv").write_text(
        "path\tbytes\tsha256\n" + "".join(f"{item['path']}\t{item['bytes']}\t{item['sha256']}\n" for item in checksums)
    )
    trace = args.outdir / "pipeline_info" / "execution_trace.txt"
    report = {
        "command": shlex.join(command),
        "pipeline": str(args.pipeline),
        "manifest": str(args.manifest),
        "reference": str(args.reference),
        "outdir": str(args.outdir),
        "host": {"system": platform.system(), "machine": platform.machine()},
        "started_at": started_at.isoformat(),
        "finished_at": finished_at.isoformat(),
        "wall_seconds": round(time.monotonic() - start, 3),
        "returncode": completed.returncode,
        "trace": trace_summary(trace),
        "software_versions": str(args.outdir / "pipeline_info" / "software_versions.yml"),
        "output_file_count": len(checksums),
        "benchmark_checksums": str(args.outdir / "benchmark_checksums.tsv"),
    }
    (args.outdir / "benchmark.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
