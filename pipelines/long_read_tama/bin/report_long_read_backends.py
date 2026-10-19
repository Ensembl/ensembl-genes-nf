#!/usr/bin/env python3
"""Build a status-first report for one or more long-read result trees.

The reporter is deliberately conservative.  A model file is not evidence that
the producing task completed: model counts are taken from explicit status
records, and unreferenced products are emitted as ``STALE_OR_UNVERIFIED``.
The output schema is stable when a run is partial, complete, resumed, or has
products from more than one attempt.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import socket
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

STATUSES = ("SUCCESS", "CACHED", "RUNNING", "FAILED", "IGNORED", "SKIPPED", "EMPTY_OUTPUT", "STALE_OR_UNVERIFIED")
BACKENDS = ("alignment", "tama", "stringtie2", "stringtie3", "tmerge", "isoquant", "bambu", "flair", "diamond", "shared")


def na(value: object = None) -> str:
    return "NA" if value is None or str(value).strip() == "" else str(value)


def normal_status(value: str) -> str:
    text = value.strip().upper().replace("-", "_")
    if text in STATUSES:
        return text
    if text in {"COMPLETE", "COMPLETED", "OK", "PASS", "VALID", "VALIDATED"}:
        return "SUCCESS"
    if text in {"CACHE", "CACHED", "CACHEABLE"}:
        return "CACHED"
    if text in {"FAIL", "ERROR", "TAMA_FAILED", "FAILED_PARTIAL"}:
        return "FAILED"
    if "IGNORE" in text:
        return "IGNORED"
    if "SKIP" in text:
        return "SKIPPED"
    if "EMPTY" in text or "NO_MODEL" in text:
        return "EMPTY_OUTPUT"
    return "STALE_OR_UNVERIFIED"


def backend_for(text: str) -> str:
    lowered = text.lower()
    for backend in BACKENDS[:-1]:
        if backend in lowered or (backend == "alignment" and any(x in lowered for x in ("minimap", "samtools", "bam_"))):
            return backend
    return "shared"


def read_table(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            delimiter = "\t" if "\t" in sample else ","
            return list(csv.DictReader(handle, delimiter=delimiter))
    except (OSError, UnicodeDecodeError):
        return []


def trace_files(root: Path) -> list[Path]:
    return [p for p in root.rglob("*") if p.is_file() and ("trace" in p.name.lower() or p.name in {"execution_trace.txt", "trace.txt"})]


def trace_rows(root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in trace_files(root):
        for row in read_table(path):
            row["_trace_path"] = str(path)
            rows.append(row)
    return rows


def status_files(root: Path, task_hashes: set[str] | None = None) -> list[Path]:
    files = []
    search_paths = [root]
    if task_hashes and root.name == "work":
        search_paths = []
        for task_hash in task_hashes:
            parts = task_hash.split("/", 1)
            if len(parts) != 2:
                continue
            search_paths.extend((root / parts[0]).glob(parts[1] + "*"))
    candidates = (path for search_path in search_paths for path in search_path.rglob("*.tsv"))
    for path in candidates:
        if not path.is_file() or "status" not in path.name.lower():
            continue
        try:
            header = path.read_text(errors="replace").splitlines()[0].lower().replace(",", "\t").split("\t")
        except (OSError, IndexError):
            continue
        # Key/value provenance manifests and model manifests also contain the
        # word "manifest", but are not task-status evidence.
        if not {"status", "state", "overall_status"}.intersection(header):
            continue
        files.append(path)
    return files


def status_rows(root: Path, task_hashes: set[str] | None = None) -> list[dict[str, str]]:
    result = []
    for path in status_files(root, task_hashes):
        for row in read_table(path):
            if not row:
                continue
            row = {str(k).strip().lower(): na(v) for k, v in row.items() if k is not None}
            row["_status_path"] = str(path)
            if task_hashes and root.name == "work":
                relative = path.relative_to(root)
                row["_trace_hash"] = "/".join(relative.parts[:2])
            raw = row.get("status", row.get("state", ""))
            if raw:
                row["status"] = normal_status(raw)
            result.append(row)
    return result


def field(row: dict[str, str], *names: str) -> str:
    for name in names:
        if row.get(name, "") not in ("", "NA", "NONE"):
            return row[name]
    return "NA"


def parse_seconds(value: str) -> float | None:
    if not value or value == "NA":
        return None
    total = 0.0
    matches = re.findall(r"([0-9]+(?:\.[0-9]+)?)\s*([dhms])", value.lower())
    for number, unit in matches:
        total += float(number) * {"d": 86400, "h": 3600, "m": 60, "s": 1}[unit]
    return total or None


def parse_memory(value: str) -> float | None:
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*([kmgt]?b)", value.lower()) if value else None
    if not match:
        return None
    return float(match.group(1)) * {"b": 1 / 1024**3, "kb": 1 / 1024**2, "mb": 1 / 1024, "gb": 1, "tb": 1024}[match.group(2)]


def repo_revision(root: Path, override: str) -> str:
    if override:
        return override
    for candidate in root.rglob("*.revision"):
        text = candidate.read_text(errors="replace").strip()
        if text:
            return text.splitlines()[0]
    return "NA"


def run_id(root: Path, override: str) -> str:
    return override or root.name


def base_row(run: dict[str, str], status: str = "STALE_OR_UNVERIFIED") -> dict[str, str]:
    return {"run_id": run["run_id"], "repo_revision": run["repo_revision"], "result_root": run["result_root"], "status": status}


def evidence_status(row: dict[str, str], traces: list[dict[str, str]]) -> str:
    status = row.get("status", "STALE_OR_UNVERIFIED")
    if status not in {"SUCCESS", "CACHED"}:
        return status
    haystack = " ".join(row.get(k, "") for k in ("backend", "process", "task_process", "_status_path"))
    candidates = [t for t in traces if row.get("_trace_hash", "") and t.get("hash", "") == row.get("_trace_hash")]
    if not candidates:
        candidates = [t for t in traces if backend_for(t.get("name", "") + " " + t.get("process", "")) == backend_for(haystack)]
    if not candidates:
        return "STALE_OR_UNVERIFIED"
    if any(str(t.get("status", "")).upper() in {"CACHED", "COMPLETED", "SUCCESS", "OK"} for t in candidates):
        return status
    return "STALE_OR_UNVERIFIED"


def discover_linked_product(row: dict[str, str], root: Path, search_roots: list[Path] | None = None) -> Path | None:
    raw = field(row, "output_path", "model_path", "bed12", "native_model", "output_file", "path")
    if raw != "NA":
        path = Path(raw)
        if not path.is_absolute():
            path = root / path
        if path.is_file():
            return path
        # Do not search a whole HPC work tree by basename: shared work trees
        # contain products from many runs.  Exact task-local discovery below is
        # the only safe fallback.
    status_path = Path(row.get("_status_path", ""))
    directory = status_path.parent if status_path else root
    candidates = [p for p in directory.iterdir()] if directory.is_dir() else []
    for path in candidates:
        if path.suffix.lower() in {".bed", ".gtf", ".gz"} and path.is_file() and path.stat().st_size:
            return path
    return None


def count_models(path: Path) -> tuple[int, int, int]:
    model_ids: set[str] = set()
    exons = 0
    junctions: set[tuple[str, int, int, str]] = set()
    try:
        opener = path.open
        lines = opener("r", errors="replace")
        for line in lines:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 12:
                chrom, start, end, name, _, strand = fields[:6]
                blocks = int(fields[9]); sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
                starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
                if blocks == len(sizes) == len(starts):
                    ex = [(int(start) + off, int(start) + off + size) for off, size in zip(starts, sizes)]
                    model_ids.add(name); exons += len(ex)
                    junctions.update((chrom, a[1], b[0], strand) for a, b in zip(ex, ex[1:]))
            elif len(fields) == 9 and fields[2] == "exon":
                attrs = fields[8]
                match = re.search(r'(?:transcript_id|transcript)\s+["=]?([^"; ]+)', attrs)
                if match:
                    model_ids.add(match.group(1))
                exons += 1
        lines.close()
    except (OSError, ValueError):
        return 0, 0, 0
    return len(model_ids), exons, len(junctions)


def bed_keys(path: Path) -> set[tuple[str, str, tuple[tuple[int, int], ...]]]:
    keys: set[tuple[str, str, tuple[tuple[int, int], ...]]] = set()
    try:
        for line in path.read_text(errors="replace").splitlines():
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.split("\t")
            if len(fields) != 12:
                continue
            chrom, start, _end, _name, _score, strand = fields[:6]
            count = int(fields[9])
            sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            if count != len(sizes) or count != len(starts):
                continue
            start_i = int(start)
            exons = tuple((start_i + off, start_i + off + size) for off, size in zip(starts, sizes))
            keys.add((chrom, strand, exons))
    except (OSError, ValueError):
        return set()
    return keys


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows([{key: na(row.get(key)) for key in fields} for row in rows])


def related_evidence_roots(root: Path) -> list[Path]:
    """Find the conventional sibling work directory for a published run."""
    names = [root.name + "-work"]
    if root.name.endswith("-results"):
        names.append(root.name[:-8] + "-work")
    if root.name == "results":
        names.append("work")
    return [candidate for name in names if (candidate := root.parent / name).is_dir() and candidate != root]


def build(root: Path, output: Path, override_run: str, override_revision: str, extra_evidence_roots: list[Path] | None = None) -> dict[str, object]:
    run = {"run_id": run_id(root, override_run), "repo_revision": repo_revision(root, override_revision), "result_root": str(root.resolve())}
    evidence_roots = [root, *related_evidence_roots(root), *(extra_evidence_roots or [])]
    traces = trace_rows(root)
    traces = [trace for trace in traces if "bos_indicus" not in (trace.get("name", "") + " " + trace.get("process", ""))]
    task_hashes = {trace.get("hash", "") for trace in traces if trace.get("hash", "") not in {"", "NA"}}
    statuses = [row for evidence_root in evidence_roots for row in status_rows(evidence_root, task_hashes)]
    statuses = [row for row in statuses if "bos_indicus" not in " ".join(row.values())]
    common = base_row(run)
    manifest = [{**common, "host": socket.gethostname(), "nextflow_version": "NA", "run_name": run["run_id"], "launch_time": "NA", "completion_time": "NA", "profile": "NA", "manifest": "NA", "reference_fasta": "NA", "minimap2_index": "NA", "parameters": "NA", "completeness": "PARTIAL_OR_UNVERIFIED"}]
    status_out: list[dict[str, str]] = []
    shard_out: list[dict[str, str]] = []
    counts: list[dict[str, str]] = []
    failures: list[dict[str, str]] = []
    merge: list[dict[str, str]] = []
    linked_products: set[Path] = set()
    for row in statuses:
        backend = field(row, "backend")
        if backend == "NA":
            trace = next((item for item in traces if item.get("hash") == row.get("_trace_hash")), {})
            backend = backend_for(field(trace, "name", "process") or field(row, "process", "_status_path"))
        scope = field(row, "scope", "level")
        unit = field(row, "run_accession", "accession", "accession_or_cohort", "cohort")
        shard = field(row, "shard", "contig")
        status = evidence_status(row, root and traces)
        product = discover_linked_product(row, root, evidence_roots)
        if product:
            linked_products.add(product.resolve())
        output_path = str(product.resolve()) if product else field(row, "output_path", "model_path", "output_file")
        base = {**common, "backend": backend, "scope": scope, "accession_or_cohort": unit, "shard": shard, "status": status, "output_path": output_path, "process": field(row, "process", "task_process"), "task_id": field(row, "task_id", "hash"), "attempt": field(row, "attempt"), "runtime": "NA", "cpu": "NA", "memory_gb": "NA", "error_or_ignore_reason": field(row, "reason", "error", "detail", "exit_status")}
        status_out.append(base)
        if shard != "NA" or field(row, "run_accession", "accession") != "NA":
            shard_out.append({**base, "mapped_reads": field(row, "mapped_reads"), "resource_class": field(row, "resource_class")})
        if status in {"FAILED", "IGNORED", "EMPTY_OUTPUT", "SKIPPED", "STALE_OR_UNVERIFIED"}:
            failures.append({**common, "backend": backend, "scope": scope, "accession_or_cohort": unit, "shard": shard, "status": status, "product_path": output_path, "evidence_path": row.get("_status_path", "NA"), "reason": base["error_or_ignore_reason"]})
        if "merge" in row.get("_status_path", "").lower() or "merge" in field(row, "process").lower():
            models = field(row, "models", "model_count", "transcript_model_count")
            merge.append({**common, "backend": backend, "scope": scope, "accession_or_cohort": unit, "status": status, "input_count": field(row, "inputs", "input_count"), "output_count": models, "input_path": field(row, "input", "input_path"), "output_path": output_path, "report_path": row.get("_status_path", "NA")})
        if product and product.suffix.lower() in {".bed", ".gtf"}:
            declared = field(row, "models", "model_count", "transcript_model_count")
            # Native validators already emit authoritative counts.  Avoid
            # rereading large per-shard GTF/BED files unless the status record
            # does not provide one.
            if declared != "NA":
                model_count, exon_count, junction_count = int(declared), "NA", "NA"
            else:
                model_count, exon_count, junction_count = "NA", "NA", "NA"
            counts.append({**common, "backend": backend, "scope": scope, "accession_or_cohort": unit, "shard": shard, "level": "shard" if shard != "NA" else ("cohort" if unit.lower() == "cohort" else "accession"), "status": status, "native_format": "bed12" if product.suffix.lower() == ".bed" else "gtf", "native_path": output_path, "normalised_bed12_path": output_path if product.suffix.lower() == ".bed" else "NA", "model_count": declared if declared != "NA" else str(model_count), "exon_count": str(exon_count), "splice_junction_count": str(junction_count), "read_assignment_path": field(row, "read_assignment_path", "read_to_model"), "tool_version": field(row, "tool_version", "version"), "container": field(row, "container", "image")})
    # A published model file may survive a resumed or superseded run.  Keep it
    # visible, but never promote it to a reportable product without evidence.
    for product in root.rglob("*"):
        if product.is_file() and product.suffix.lower() in {".bed", ".gtf"} and product.resolve() not in linked_products:
            backend = backend_for(str(product))
            stale = {**common, "backend": backend, "scope": "NA", "accession_or_cohort": "NA", "shard": "NA", "status": "STALE_OR_UNVERIFIED", "output_path": str(product.resolve()), "process": "NA", "task_id": "NA", "attempt": "NA", "runtime": "NA", "cpu": "NA", "memory_gb": "NA", "error_or_ignore_reason": "model file has no linked status record"}
            status_out.append(stale)
            failures.append({**common, "backend": backend, "scope": "NA", "accession_or_cohort": "NA", "shard": "NA", "status": stale["status"], "product_path": stale["output_path"], "evidence_path": "NA", "reason": stale["error_or_ignore_reason"]})
    resources: list[dict[str, str]] = []
    for trace in traces:
        name = field(trace, "name", "process")
        status = normal_status(field(trace, "status"))
        if status == "STALE_OR_UNVERIFIED" and field(trace, "status") == "NA":
            status = "RUNNING"
        runtime = parse_seconds(field(trace, "realtime", "duration")); memory = parse_memory(field(trace, "peak_rss", "peak_rss_gb"))
        resources.append({**common, "backend": backend_for(name), "process": name, "status": status, "task_id": field(trace, "task_id", "hash"), "attempt": field(trace, "attempt"), "runtime_seconds": str(runtime) if runtime is not None else "NA", "cpu_hours": field(trace, "cpus", "cpu") if runtime is None else str((runtime * float(field(trace, "cpus", "cpu"))) / 3600) if field(trace, "cpus", "cpu") != "NA" else "NA", "peak_rss_gb": str(memory) if memory is not None else "NA", "cache_status": "CACHED" if "cache" in field(trace, "status").lower() else "NA", "trace_path": trace.get("_trace_path", "NA")})
    for row in status_out:
        if row["status"] in {"SUCCESS", "CACHED"} and row["output_path"] == "NA":
            row["status"] = "STALE_OR_UNVERIFIED"
    run_dir = output / run["run_id"]
    fields_common = ["run_id", "repo_revision", "result_root", "status"]
    write_tsv(run_dir / "run_manifest.tsv", fields_common + ["host", "nextflow_version", "run_name", "launch_time", "completion_time", "profile", "manifest", "reference_fasta", "minimap2_index", "parameters", "completeness"], manifest)
    write_tsv(run_dir / "input_alignment_summary.tsv", fields_common + ["accession", "read_count", "mapped_reads", "primary_mapped_reads", "supplementary_reads", "mapping_rate", "reference", "shard_count", "skipped_shard_keys"], [])
    write_tsv(run_dir / "backend_status.tsv", fields_common + ["backend", "scope", "accession_or_cohort", "shard", "output_path", "process", "task_id", "attempt", "runtime", "cpu", "memory_gb", "error_or_ignore_reason"], status_out)
    write_tsv(run_dir / "shard_status.tsv", fields_common + ["backend", "scope", "accession_or_cohort", "shard", "output_path", "process", "task_id", "attempt", "runtime", "cpu", "memory_gb", "error_or_ignore_reason", "mapped_reads", "resource_class"], shard_out)
    write_tsv(run_dir / "backend_model_counts.tsv", fields_common + ["backend", "scope", "accession_or_cohort", "shard", "level", "native_format", "native_path", "normalised_bed12_path", "model_count", "exon_count", "splice_junction_count", "read_assignment_path", "tool_version", "container"], counts)
    comparisons: list[dict[str, str]] = []
    comparable: dict[tuple[str, str, str], dict[str, set]] = {}
    for row in counts:
        if row["status"] not in {"SUCCESS", "CACHED"} or row["native_format"] != "bed12" or row["native_path"] == "NA":
            continue
        path = Path(row["native_path"])
        comparable.setdefault((row["scope"], row["accession_or_cohort"], row["level"]), {})[row["backend"]] = bed_keys(path)
    for (scope, unit, level), backend_sets in comparable.items():
        names = sorted(backend_sets)
        for index, left in enumerate(names):
            for right in names[index + 1:]:
                left_set, right_set = backend_sets[left], backend_sets[right]
                union = left_set | right_set
                shared = left_set & right_set
                comparisons.append({**common, "status": "SUCCESS", "left_backend": left, "right_backend": right, "scope": scope, "accession_or_cohort": unit, "shared_models": str(len(shared)), "left_specific_models": str(len(left_set - right_set)), "right_specific_models": str(len(right_set - left_set)), "jaccard": str(len(shared) / len(union)) if union else "NA", "denominator": str(len(union)), "filtering": f"exact BED12 exon-chain key; level={level}; only SUCCESS/CACHED products"})
    write_tsv(run_dir / "comparison_summary.tsv", fields_common + ["left_backend", "right_backend", "scope", "accession_or_cohort", "shared_models", "left_specific_models", "right_specific_models", "jaccard", "denominator", "filtering"], comparisons)
    write_tsv(run_dir / "merge_summary.tsv", fields_common + ["backend", "scope", "accession_or_cohort", "input_count", "output_count", "input_path", "output_path", "report_path"], merge)
    write_tsv(run_dir / "failure_summary.tsv", fields_common + ["backend", "scope", "accession_or_cohort", "shard", "product_path", "evidence_path", "reason"], failures)
    write_tsv(run_dir / "resource_summary.tsv", fields_common + ["backend", "process", "task_id", "attempt", "runtime_seconds", "cpu_hours", "peak_rss_gb", "cache_status", "trace_path"], resources)
    write_tsv(run_dir / "software_versions.tsv", fields_common + ["backend", "process", "tool_version", "container", "evidence_path"], [])
    successful = sum(row["status"] in {"SUCCESS", "CACHED"} for row in status_out)
    md = [f"# Long-read backend report: {run['run_id']}", "", f"Result root: `{run['result_root']}`", f"Repository revision: `{run['repo_revision']}`", "", "This report is status-first. Products without a matching status record and trace evidence are labelled `STALE_OR_UNVERIFIED`; unknown values remain `NA`.", "", f"Execution/status records: {len(status_out)}; reportable success/cached records: {successful}; trace rows: {len(traces)}.", "", "## Backend status", "", "| Backend | Scope | Unit | Shard | Status | Output |", "|---|---|---|---|---|---|"]
    md.extend(f"| {r['backend']} | {r['scope']} | {r['accession_or_cohort']} | {r['shard']} | {r['status']} | `{r['output_path']}` |" for r in status_out[:200])
    md += ["", "## Known interpretation cautions", "", "- FLAIR products are provisional unless their status evidence identifies the corrected Wave container and contig-preserving junction extraction.", "- Bambu cohort products from the affected parallel reducer attempt remain failed/ignored until a successful single-core reducer status is present.", "- TAMA counts are provisional while any expected collapse shard is running, skipped, failed, or unverified.", "- A `skip_tama` directory name is not used as evidence of skipping; statuses and traces are authoritative.", "", "## Machine-readable tables", ""]
    md.extend(f"- [{name}](./{name})" for name in ("run_manifest.tsv", "input_alignment_summary.tsv", "backend_status.tsv", "shard_status.tsv", "backend_model_counts.tsv", "comparison_summary.tsv", "merge_summary.tsv", "failure_summary.tsv", "resource_summary.tsv", "software_versions.tsv"))
    (run_dir / "report.md").write_text("\n".join(md) + "\n")
    return {"run_id": run["run_id"], "output": str(run_dir), "status_rows": len(status_out), "trace_rows": len(traces), "model_rows": len(counts)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", action="append", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--run-id")
    parser.add_argument("--repo-revision", default="")
    parser.add_argument("--evidence-root", action="append", type=Path, default=[], help="Additional work/status tree used to link published products to trace evidence")
    args = parser.parse_args()
    summaries = []
    for root in args.results_root:
        if not root.is_dir():
            raise SystemExit(f"results root is not a directory: {root}")
        summaries.append(build(root, args.output_dir, args.run_id or "", args.repo_revision, args.evidence_root))
    (args.output_dir / "index.json").write_text(json.dumps({"reports": summaries}, indent=2) + "\n")


if __name__ == "__main__":
    main()
