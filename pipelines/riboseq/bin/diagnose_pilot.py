#!/usr/bin/env python3
"""Explain per-run outcomes from a published Ribo-seq pilot.

This is intentionally a read-only report generator.  It joins the RDP
selection table with published getRPF reports, QC metrics, gate decisions,
and sample completeness.  A run can therefore be distinguished between
structure inference failure, missing/software output, and genuine QC-gate
failure instead of treating every absent downstream file as bad biology.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


TRUE = {"1", "true", "yes", "y", "pass", "passed", "ok"}


def as_bool(value):
    return str(value or "").strip().lower() in TRUE


def read_tsv(path):
    if not path or not path.exists() or not path.stat().st_size:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_selection(path):
    rows = read_tsv(path)
    result = {}
    for row in rows:
        run = row.get("run_id") or row.get("Run") or row.get("run") or row.get("accession")
        if not run:
            continue
        result[run] = {
            "run_id": run,
            "study_accession": row.get("study_accession") or row.get("study") or "",
            "organism": row.get("organism") or row.get("species") or "",
            "depth_bin": row.get("depth_bin") or "",
        }
    return result


def sample_from_name(path):
    name = path.name
    sample = None
    for suffix in ("_extraction_report.json", ".extraction_report.json", "_qc_metrics.tsv", "_qc_artifacts.tsv", "_gate_selection.tsv", "_qc_eval.tsv"):
        if name.endswith(suffix):
            sample = name[: -len(suffix)]
            break
    if sample is None and "_RiboMetric" in name:
        sample = name.split("_RiboMetric", 1)[0]
    elif sample is None:
        sample = path.stem.split(".", 1)[0]
    # COLLECT_QC_METRICS and QC_GATE include the Nextflow task index in the
    # published filename (for example SRR123.1.qc_metrics.tsv).
    sample = re.sub(r"_trimmed$", "", sample)
    return re.sub(r"\.\d+$", "", sample)


def load_extraction(root):
    reports = {}
    for path in sorted((root / "getRPF" / "extract").glob("*extraction_report.json")):
        sample = sample_from_name(path)
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            payload = {}
        report = payload.get("results", payload) if isinstance(payload, dict) else {}
        if not isinstance(report, dict):
            report = {}
        reasons = report.get("withheld_reasons") or report.get("withhold_reasons") or report.get("reasons") or []
        if isinstance(reasons, str):
            reasons = [reasons]
        status = str(report.get("status") or report.get("decision") or "").lower()
        emitted = report.get("emitted")
        if emitted is True or status in {"emitted", "pass", "passed", "ok"}:
            outcome = "emitted"
        elif emitted is False or status in {"withheld", "failed", "fail", "rejected"} or reasons:
            outcome = "withheld"
        else:
            outcome = "unknown"
        reports[sample] = {
            "extraction_outcome": outcome,
            "extraction_reasons": ";".join(str(item) for item in reasons),
            "extraction_path": str(path),
        }
    return reports


def load_metrics(root):
    metrics = {}
    for path in sorted((root / "qc_gate").glob("*_qc_metrics.tsv")):
        sample = sample_from_name(path)
        row = metrics.setdefault(sample, {"qc_metrics_path": str(path)})
        for metric in read_tsv(path):
            name = metric.get("metric", "")
            value = metric.get("value_num", "")
            if name and value != "":
                row[name] = value
    return metrics


def load_gate(root):
    gate = {}
    for path in sorted((root / "qc_gate").glob("*_gate_selection.tsv")):
        sample = sample_from_name(path)
        row = gate.setdefault(sample, {"gate_path": str(path), "gate_selected": False, "gate_tier": ""})
        for decision in read_tsv(path):
            if decision.get("tier") == "good":
                row["gate_selected"] = as_bool(decision.get("selected_for_translon")) or as_bool(decision.get("selected_for_trackhub"))
                row["gate_tier"] = decision.get("tier", "")
                row["gate_reason"] = decision.get("selected_reason", "")
                row["pass_lengths"] = decision.get("pass_lengths_json", "")
    return gate


def load_failed_rules(root):
    failed = {}
    for path in sorted((root / "qc_gate").glob("*_qc_eval.tsv")):
        sample = sample_from_name(path)
        names = failed.setdefault(sample, [])
        for row in read_tsv(path):
            if not as_bool(row.get("passed")) and row.get("rule_name"):
                names.append(row["rule_name"])
    return {sample: ";".join(sorted(set(names))) for sample, names in failed.items()}


def load_status(root):
    paths = [root / "pipeline_info" / "sample_status.tsv", root / "sample_status.tsv"]
    for path in paths:
        rows = read_tsv(path)
        if rows:
            return {
                row.get("sample_id") or row.get("run_id") or row.get("Run"): row
                for row in rows
                if row.get("sample_id") or row.get("run_id") or row.get("Run")
            }
    return {}


def classify(row):
    if row.get("extraction_outcome") == "withheld":
        return "withheld_at_structure_inference"
    if row.get("gate_selected"):
        return "pass_good_gate"
    if row.get("extraction_outcome") == "emitted" and not row.get("qc_metrics_path"):
        return "emitted_but_qc_missing_or_software_failure"
    if row.get("extraction_outcome") == "emitted" and row.get("qc_metrics_path"):
        return "emitted_but_qc_gate_failed"
    if row.get("status"):
        return "present_in_status_only"
    return "unresolved"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-root", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    selection = read_selection(args.selection)
    extraction = load_extraction(args.pilot_root)
    metrics = load_metrics(args.pilot_root)
    gate = load_gate(args.pilot_root)
    failed_rules = load_failed_rules(args.pilot_root)
    status = load_status(args.pilot_root)

    samples = set(selection) | set(extraction) | set(metrics) | set(gate) | set(status)
    output_rows = []
    for sample in sorted(samples):
        row = dict(selection.get(sample, {"run_id": sample}))
        row.update(extraction.get(sample, {}))
        row.update(metrics.get(sample, {}))
        row.update(gate.get(sample, {}))
        row["failed_qc_rules"] = failed_rules.get(sample, "")
        row.update({f"status_{key}": value for key, value in status.get(sample, {}).items()})
        row["status"] = status.get(sample, {}).get("status", "")
        row["diagnosis"] = classify(row)
        output_rows.append(row)

    preferred = [
        "run_id", "study_accession", "organism", "depth_bin", "diagnosis",
        "extraction_outcome", "extraction_reasons", "status",
        "gate_selected", "gate_reason", "failed_qc_rules",
        "star.unique_mapped_pct", "getrpf.checks_pass_frac",
        "ribometric.periodicity_dominance", "ribometric.information_content",
        "ribometric.n_recommended_read_lengths", "ribometric.recommended_read_proportion",
        "qc_metrics_path", "extraction_path", "gate_path",
    ]
    keys = list(preferred)
    for row in output_rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)

    counts = Counter(row["diagnosis"] for row in output_rows)
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps({"n_runs": len(output_rows), "diagnoses": counts}, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"n_runs": len(output_rows), "diagnoses": counts}, sort_keys=True))


if __name__ == "__main__":
    main()
