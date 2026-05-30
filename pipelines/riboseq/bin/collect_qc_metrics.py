#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from pathlib import Path


STAR_KEYS = [
    ("Number of input reads", "star.total_reads"),
    ("Uniquely mapped reads number", "star.unique_mapped_reads"),
    ("Uniquely mapped reads %", "star.unique_mapped_pct"),
    ("Number of reads mapped to multiple loci", "star.multi_mapped_reads"),
    ("% of reads mapped to multiple loci", "star.multi_mapped_pct"),
]

PASS_VALUES = {"true", "pass", "passed", "1", "ok", "yes", "y"}


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def add_artifact(rows, run_id, sample_id, module, label, path, content_type):
    path = Path(path)
    if not path.exists():
        return
    rows.append(
        {
            "artifact_id": f"{sample_id}:{module}:{label}",
            "run_id": run_id,
            "sample_id": sample_id,
            "module": module,
            "name": path.name,
            "path": str(path),
            "checksum": md5sum(path),
            "size_bytes": path.stat().st_size,
            "content_type": content_type,
        }
    )


def parse_star_log(path):
    vals = {}
    with open(path) as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = [part.strip() for part in line.split("|", 1)]
            for star_key, metric_name in STAR_KEYS:
                if key != star_key:
                    continue
                try:
                    vals[metric_name] = float(value.replace("%", "").strip())
                except ValueError:
                    pass

    total = vals.get("star.total_reads") or 0.0
    unique = vals.get("star.unique_mapped_reads") or 0.0
    if total:
        vals["star.unique_mapped_frac"] = unique / total
    return vals


def parse_getrpf_checks(path):
    ok = 0.0
    total = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            total += 1
            if parts[1].strip().lower() in PASS_VALUES:
                ok += 1.0
    return (ok / total) if total else None


def load_json(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return {}


def load_csv_rows(path):
    if not path or not Path(path).exists():
        return []
    with open(path) as handle:
        return list(csv.DictReader(handle))


def add_metric(rows, run_id, sample_id, module, metric, scope, length, value, unit=None):
    try:
        value_num = float(value)
    except (TypeError, ValueError):
        return
    rows.append(
        {
            "run_id": run_id,
            "sample_id": sample_id,
            "module": module,
            "metric": metric,
            "scope": scope,
            "length": "" if length is None else int(length),
            "value_num": value_num,
            "unit": "" if unit is None else unit,
        }
    )


def collect_star(rows, args):
    for metric, value in parse_star_log(args.star_log).items():
        add_metric(rows, args.run_id, args.sample_id, "STAR", metric, "sample", None, value)


def collect_getrpf(rows, args):
    if not args.getrpf_report or not args.getrpf_checks:
        return
    if not Path(args.getrpf_report).stat().st_size or not Path(args.getrpf_checks).stat().st_size:
        return
    pass_fraction = parse_getrpf_checks(args.getrpf_checks)
    if pass_fraction is not None:
        add_metric(
            rows,
            args.run_id,
            args.sample_id,
            "getRPF",
            "getrpf.checks_pass_frac",
            "sample",
            None,
            pass_fraction,
        )


def collect_ribometric(rows, args):
    jm = load_json(args.ribometric_json) if args.ribometric_json else {}
    ic = None
    for key in ("information_content", "info_content"):
        if key in jm:
            ic = jm[key]
            break
    if ic is None:
        try:
            ic = jm.get("metagene", {}).get("information_content")
        except Exception:
            pass
    if ic is not None:
        add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.information_content", "sample", None, ic)

    f0 = jm.get("frame0_fraction") or (jm.get("frames", {}).get("frame0") if isinstance(jm.get("frames"), dict) else None)
    if f0 is not None:
        add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.frame0_frac_overall", "sample", None, f0)

    per_len = jm.get("per_length") or jm.get("length_metrics") or {}
    if per_len:
        for key, values in per_len.items():
            try:
                read_len = int(key)
            except Exception:
                continue
            reads = values.get("reads") or values.get("count") or values.get("n")
            if reads is not None:
                add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.reads_at_length", "length", read_len, reads, "reads")
            frame0 = values.get("frame0_fraction") or values.get("frame0") or values.get("f0_frac")
            if frame0 is not None:
                add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.frame0_frac", "length", read_len, frame0)
        return

    for row in load_csv_rows(args.ribometric_csv):
        try:
            read_len = int(row.get("length") or row.get("read_length"))
        except Exception:
            continue
        total = row.get("total") or row.get("reads") or row.get("n")
        frame0 = row.get("frame0_frac") or row.get("f0") or row.get("frame0_fraction")
        if frame0 is None and "frame0" in row and total:
            try:
                frame0 = float(row["frame0"]) / float(total)
            except Exception:
                pass
        if total is not None:
            add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.reads_at_length", "length", read_len, total, "reads")
        if frame0 is not None:
            add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.frame0_frac", "length", read_len, frame0)


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--study-id", default="unknown")
    parser.add_argument("--star-log", required=True)
    parser.add_argument("--getrpf-report", required=True)
    parser.add_argument("--getrpf-checks", required=True)
    parser.add_argument("--ribometric-json", required=True)
    parser.add_argument("--ribometric-csv", required=True)
    parser.add_argument("--offsets", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    metrics = []
    artifacts = []

    collect_star(metrics, args)
    collect_getrpf(metrics, args)
    collect_ribometric(metrics, args)

    add_artifact(artifacts, args.run_id, args.sample_id, "STAR", "Log.final.out", args.star_log, "text/plain")
    add_artifact(artifacts, args.run_id, args.sample_id, "getRPF", "getRPF_report", args.getrpf_report, "text/plain")
    add_artifact(artifacts, args.run_id, args.sample_id, "getRPF", "getRPF_checks", args.getrpf_checks, "text/plain")
    add_artifact(artifacts, args.run_id, args.sample_id, "RiboMetric", "json", args.ribometric_json, "application/json")
    add_artifact(artifacts, args.run_id, args.sample_id, "RiboMetric", "csv", args.ribometric_csv, "text/csv")
    add_artifact(artifacts, args.run_id, args.sample_id, "RiboMetric", "offsets", args.offsets, "text/tab-separated-values")

    write_tsv(
        f"{args.out_prefix}.qc_metrics.tsv",
        metrics,
        ["run_id", "sample_id", "module", "metric", "scope", "length", "value_num", "unit"],
    )
    write_tsv(
        f"{args.out_prefix}.qc_artifacts.tsv",
        artifacts,
        ["artifact_id", "run_id", "sample_id", "module", "name", "path", "checksum", "size_bytes", "content_type"],
    )


if __name__ == "__main__":
    main()
