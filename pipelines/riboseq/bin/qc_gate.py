#!/usr/bin/env python3
import argparse
import csv
import json
import shutil
from pathlib import Path

import yaml


OPS = {
    ">=": lambda observed, threshold: observed is not None and observed >= threshold,
    "<=": lambda observed, threshold: observed is not None and observed <= threshold,
    "between": lambda observed, threshold: observed is not None and threshold[0] <= observed <= threshold[1],
}


def read_metrics(path):
    metrics = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            length = int(row["length"]) if row["length"] else None
            key = (row["metric"], row["scope"], length)
            try:
                metrics[key] = float(row["value_num"])
            except ValueError:
                continue
    return metrics


def get_metric(metrics, metric, scope, length=None):
    return metrics.get((metric, scope, int(length) if length is not None else None))


def read_offsets_table(path):
    rows = []
    with open(path) as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                read_len = int(fields[0])
                offset = int(fields[1])
            except (IndexError, ValueError):
                continue
            rows.append((read_len, offset))
    return rows


def write_filtered_offsets(in_path, out_path, keep_lengths):
    with open(in_path) as source, open(out_path, "w") as dest:
        header = source.readline()
        dest.write(header)
        for line in source:
            if not line.strip():
                continue
            try:
                read_len = int(line.split("\t", 1)[0])
            except ValueError:
                continue
            if read_len in keep_lengths:
                dest.write(line)


def evaluate_rule(metrics, args, rule, rule_set_id, scope, length=None):
    metric = rule["metric"]
    op = rule["op"]
    threshold = rule["threshold"]
    observed = get_metric(metrics, metric, scope, length)
    if observed is None and not rule.get("include_if_missing", False):
        passed = False
    elif isinstance(threshold, (list, tuple)) and op == "between":
        passed = OPS["between"](observed, threshold)
    else:
        passed = OPS[op](observed, threshold)

    threshold_min = float(threshold[0]) if isinstance(threshold, (list, tuple)) else float(threshold)
    threshold_max = float(threshold[1]) if isinstance(threshold, (list, tuple)) and len(threshold) > 1 else None
    return passed, (
        args.run_id,
        args.sample_id,
        rule_set_id,
        scope,
        "" if length is None else int(length),
        rule["name"],
        metric,
        op,
        threshold_min,
        "" if threshold_max is None else threshold_max,
        "" if observed is None else float(observed),
        "",
        rule["severity"],
        bool(passed),
    )


def write_qc_eval(path, records):
    fieldnames = [
        "run_id",
        "sample_id",
        "rule_set_id",
        "scope",
        "length",
        "rule_name",
        "metric",
        "op",
        "threshold_num_min",
        "threshold_num_max",
        "observed_num",
        "observed_text",
        "severity",
        "passed",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fieldnames)
        writer.writerows(records)


def write_gate_selection(path, args, rule_set_id, selected_for_translon, selected_for_trackhub, reason, keep_lengths):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "run_id",
                "sample_id",
                "rule_set_id",
                "selected_for_translon",
                "selected_for_trackhub",
                "selected_reason",
                "pass_lengths_json",
            ]
        )
        writer.writerow(
            [
                args.run_id,
                args.sample_id,
                rule_set_id,
                bool(selected_for_translon),
                bool(selected_for_trackhub),
                reason,
                json.dumps(sorted(keep_lengths)),
            ]
        )


def write_rule_set(path, rule_set_id, name, version, yaml_text):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rule_set_id", "name", "version", "yaml_text"])
        writer.writerow([rule_set_id, name, int(version), yaml_text])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--metrics-tsv", required=True)
    parser.add_argument("--best-offset", required=True)
    parser.add_argument("--rules-yaml", required=True)
    parser.add_argument("--rule-set-name", default="default")
    parser.add_argument("--apply-for", default="translon,trackhub")
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    metrics = read_metrics(args.metrics_tsv)
    rules_text = Path(args.rules_yaml).read_text()
    rules = yaml.safe_load(rules_text) or {}
    rule_set_name = rules.get("rule_set_name", args.rule_set_name)
    rule_set_version = int(rules.get("version", 1))
    rule_set_id = f"{rule_set_name}_v{rule_set_version}"

    keep_lengths = set()
    eval_records = []

    for read_len, _offset in read_offsets_table(args.best_offset):
        passed_all = True
        for rule in rules.get("length_rules", []):
            passed, record = evaluate_rule(metrics, args, rule, rule_set_id, "length", read_len)
            eval_records.append(record)
            if not passed:
                passed_all = False
        if passed_all:
            keep_lengths.add(read_len)

    sample_pass = True
    for rule in rules.get("sample_rules", []):
        passed, record = evaluate_rule(metrics, args, rule, rule_set_id, "sample")
        eval_records.append(record)
        if not passed:
            sample_pass = False

    apply_for = {item.strip() for item in args.apply_for.split(",") if item.strip()}
    selected = sample_pass and len(keep_lengths) > 0
    selected_for_translon = "translon" in apply_for and selected
    selected_for_trackhub = "trackhub" in apply_for and selected
    reason = "sample rules pass and at least one passing length" if selected else "rules not met"

    with open(f"{args.out_prefix}.pass_lengths.tsv", "w") as handle:
        handle.write("length\n")
        for read_len in sorted(keep_lengths):
            handle.write(f"{read_len}\n")

    pass_offsets = f"{args.out_prefix}.offsets.pass.tsv"
    write_filtered_offsets(args.best_offset, pass_offsets, keep_lengths)
    if selected_for_translon or selected_for_trackhub:
        shutil.copyfile(pass_offsets, f"{args.out_prefix}.offsets.selected.tsv")
    if selected_for_translon:
        Path(f"{args.out_prefix}.translon.selected.txt").write_text(f"{args.sample_id}\n")

    write_qc_eval(f"{args.out_prefix}.qc_eval.tsv", eval_records)
    write_gate_selection(
        f"{args.out_prefix}.gate_selection.tsv",
        args,
        rule_set_id,
        selected_for_translon,
        selected_for_trackhub,
        reason,
        keep_lengths,
    )
    write_rule_set(f"{args.out_prefix}.qc_rule_set.tsv", rule_set_id, rule_set_name, rule_set_version, rules_text)

    with open(f"{args.out_prefix}.qc.json", "w") as handle:
        json.dump(
            {
                "sample_id": args.sample_id,
                "rule_set_id": rule_set_id,
                "selected_for_translon": selected_for_translon,
                "selected_for_trackhub": selected_for_trackhub,
                "selected_reason": reason,
                "pass_lengths": sorted(keep_lengths),
            },
            handle,
        )


if __name__ == "__main__":
    main()
