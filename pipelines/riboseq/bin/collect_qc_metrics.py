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
    with open(path, newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            dialect = csv.excel
        return list(csv.DictReader(handle, dialect=dialect))


def results_root(report):
    return report.get("results", report) if isinstance(report, dict) else {}


def nested_metric(metrics, key, nested_key=None):
    value = metrics.get(key)
    if nested_key is not None and isinstance(value, dict):
        return value.get(nested_key)
    return value


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


def has_metric(rows, metric, scope, length=None):
    length_value = "" if length is None else int(length)
    return any(
        row["metric"] == metric and row["scope"] == scope and row["length"] == length_value
        for row in rows
    )


def add_metric_if_missing(rows, run_id, sample_id, module, metric, scope, length, value, unit=None):
    if not has_metric(rows, metric, scope, length):
        add_metric(rows, run_id, sample_id, module, metric, scope, length, value, unit)


def parse_read_len(row):
    try:
        value = row.get("length") or row.get("read_length") or row.get("read_len")
        return int(float(value))
    except Exception:
        return None


def parse_ribometric_metric_name(metric_name):
    if not metric_name:
        return None, None
    for suffix in ("_global", "_sample"):
        if metric_name.endswith(suffix):
            return metric_name[: -len(suffix)], None
    for separator in ("_rl", "_"):
        prefix, sep, suffix = metric_name.rpartition(separator)
        if not sep:
            continue
        try:
            return prefix, int(float(suffix))
        except ValueError:
            continue
    return metric_name, None


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
    report = results_root(load_json(args.ribometric_json)) if args.ribometric_json else {}
    metrics = report.get("metrics", {}) if isinstance(report.get("metrics"), dict) else {}

    ic = None
    for key in ("information_content", "info_content", "periodicity_information_weighted_score"):
        if key in report:
            ic = report[key]
            break
        if key in metrics and not isinstance(metrics[key], dict):
            ic = metrics[key]
            break
    if ic is None:
        ic = nested_metric(metrics, "periodicity_information", "global")
    if ic is None and isinstance(report.get("metagene"), dict):
        ic = report["metagene"].get("information_content")
    if ic is not None:
        add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.information_content", "sample", None, ic)

    dominance = metrics.get("periodicity_dominance")
    if isinstance(dominance, dict):
        for key, value in dominance.items():
            if str(key) == "global":
                add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.periodicity_dominance", "sample", None, value)
                continue
            try:
                read_len = int(float(key))
            except Exception:
                continue
            add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.periodicity_dominance", "length", read_len, value)

    recommended = report.get("recommended_read_lengths")
    n_recommended = (
        report.get("n_recommended_read_lengths")
        if report.get("n_recommended_read_lengths") is not None
        else metrics.get("n_recommended_read_lengths")
    )
    recommended_read_proportion = (
        report.get("recommended_read_proportion")
        if report.get("recommended_read_proportion") is not None
        else metrics.get("recommended_read_proportion")
    )

    if not isinstance(recommended, dict):
        recommended = derive_recommended_read_lengths(report)
        if n_recommended is None:
            n_recommended = recommended.get("n_recommended")
        if recommended_read_proportion is None:
            recommended_read_proportion = recommended.get("recommended_read_proportion")

    if n_recommended is not None:
        add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.n_recommended_read_lengths", "sample", None, n_recommended)
    if recommended_read_proportion is not None:
        add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.recommended_read_proportion", "sample", None, recommended_read_proportion)
    collect_recommended_read_length_metrics(rows, args, recommended)

    f0 = report.get("frame0_fraction") or (report.get("frames", {}).get("frame0") if isinstance(report.get("frames"), dict) else None)
    if f0 is not None:
        add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.frame0_frac_overall", "sample", None, f0)

    read_lengths = report.get("read_length_distribution", {})
    read_frames = report.get("read_frame_distribution", {})
    if isinstance(read_lengths, dict):
        for key, count in read_lengths.items():
            try:
                read_len = int(float(key))
            except Exception:
                continue
            add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.reads_at_length", "length", read_len, count, "reads")

    if isinstance(read_frames, dict):
        for key, frames in read_frames.items():
            try:
                read_len = int(float(key))
            except Exception:
                continue
            if not isinstance(frames, dict):
                continue
            try:
                frame0 = float(frames.get("0", frames.get(0, 0)))
                total = sum(float(frames.get(str(frame), frames.get(frame, 0))) for frame in (0, 1, 2))
            except (TypeError, ValueError):
                continue
            if total:
                add_metric(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.frame0_frac", "length", read_len, frame0 / total)

    if read_lengths or read_frames:
        collect_ribometric_csv_supplements(rows, args)
        return

    per_len = report.get("per_length") or report.get("length_metrics") or {}
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
        collect_ribometric_csv_supplements(rows, args)
        return

    for row in load_csv_rows(args.ribometric_csv):
        try:
            read_len = parse_read_len(row)
            if read_len is None:
                continue
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
            add_metric_if_missing(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.reads_at_length", "length", read_len, total, "reads")
        if frame0 is not None:
            add_metric_if_missing(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.frame0_frac", "length", read_len, frame0)

    collect_ribometric_csv_supplements(rows, args)


def collect_ribometric_csv_supplements(rows, args):
    for row in iter_metric_score_rows(args.ribometric_csv):
        metric = row.get("metric")
        score = row.get("score")
        metric_base, read_len = parse_ribometric_metric_name(metric)
        if metric_base in {"periodicity_information", "periodicity_information_weighted_score"}:
            add_metric_if_missing(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.information_content", "sample", None, score)
        elif metric_base == "periodicity_dominance":
            scope = "sample" if read_len is None else "length"
            add_metric_if_missing(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.periodicity_dominance", scope, read_len, score)
        elif metric_base == "recommended_read_proportion":
            add_metric_if_missing(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.recommended_read_proportion", "sample", None, score)
        elif metric_base == "n_recommended_read_lengths":
            add_metric_if_missing(rows, args.run_id, args.sample_id, "RiboMetric", "ribometric.n_recommended_read_lengths", "sample", None, score)


def iter_metric_score_rows(path):
    if not path or not Path(path).exists():
        return []
    rows = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        all_rows = [row for row in reader if row]
    if not all_rows:
        return rows

    first = [cell.strip().lower() for cell in all_rows[0]]
    if "metric" in first:
        metric_idx = first.index("metric")
        score_idx = (
            first.index("score")
            if "score" in first
            else first.index("value")
            if "value" in first
            else 1
        )
        data_rows = all_rows[1:]
    else:
        metric_idx = 0
        score_idx = 1
        data_rows = all_rows

    for row in data_rows:
        if len(row) <= max(metric_idx, score_idx):
            continue
        rows.append({"metric": row[metric_idx], "score": row[score_idx]})
    return rows


def derive_recommended_read_lengths(report, min_periodicity=0.5, min_read_proportion=0.05):
    read_lengths = report.get("read_length_distribution", {})
    read_frames = report.get("read_frame_distribution", {})
    if not isinstance(read_lengths, dict) or not isinstance(read_frames, dict):
        return {}

    try:
        total_reads = sum(float(count) for count in read_lengths.values())
    except (TypeError, ValueError):
        return {}
    if not total_reads:
        return {}

    by_read_length = {}
    recommended_lengths = []
    for key, frames in read_frames.items():
        if not isinstance(frames, dict):
            continue
        try:
            read_len = int(float(key))
            frame_counts = [float(frames.get(str(frame), frames.get(frame, 0))) for frame in (0, 1, 2)]
            frame_total = sum(frame_counts)
            read_count = float(read_lengths.get(str(read_len), read_lengths.get(read_len, 0)))
        except (TypeError, ValueError):
            continue
        periodicity = max(frame_counts) / frame_total if frame_total else 0.0
        read_proportion = read_count / total_reads
        is_recommended = periodicity >= min_periodicity and read_proportion >= min_read_proportion
        by_read_length[str(read_len)] = {
            "periodicity": periodicity,
            "read_proportion": read_proportion,
            "recommended": is_recommended,
        }
        if is_recommended:
            recommended_lengths.append(read_len)

    recommended_read_proportion = sum(
        float(read_lengths.get(str(read_len), read_lengths.get(read_len, 0)))
        for read_len in recommended_lengths
    ) / total_reads
    return {
        "by_read_length": by_read_length,
        "recommended_lengths": sorted(recommended_lengths),
        "n_recommended": len(recommended_lengths),
        "recommended_read_proportion": recommended_read_proportion,
    }


def collect_recommended_read_length_metrics(rows, args, recommended):
    if not isinstance(recommended, dict):
        return
    by_read_length = recommended.get("by_read_length", {})
    if not isinstance(by_read_length, dict):
        return
    for key, values in by_read_length.items():
        if not isinstance(values, dict):
            continue
        try:
            read_len = int(float(key))
        except (TypeError, ValueError):
            continue
        if values.get("recommended") is not None:
            add_metric_if_missing(
                rows,
                args.run_id,
                args.sample_id,
                "RiboMetric",
                "ribometric.recommended_length",
                "length",
                read_len,
                1.0 if values.get("recommended") else 0.0,
            )
        if values.get("read_proportion") is not None:
            add_metric_if_missing(
                rows,
                args.run_id,
                args.sample_id,
                "RiboMetric",
                "ribometric.read_proportion",
                "length",
                read_len,
                values.get("read_proportion"),
            )


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
