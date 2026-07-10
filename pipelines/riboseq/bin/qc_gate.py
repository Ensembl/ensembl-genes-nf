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


def collect_named_rules(rules):
    named = {}
    for section in ("sample_rules", "length_rules", "additional_rules"):
        for rule in rules.get(section, []):
            named[rule["name"]] = rule
    return named


def resolve_tier_rules(rules, tier_config, scope):
    defaults = rules.get(f"{scope}_rules", [])
    requested = (tier_config or {}).get(f"{scope}_rules")
    if not requested:
        return defaults
    named = collect_named_rules(rules)
    return [named[name] for name in requested if name in named]


def evaluate_sample_rules(metrics, args, sample_rules, rule_set_id, eval_records):
    passed_all = True
    for rule in sample_rules:
        passed, record = evaluate_rule(metrics, args, rule, rule_set_id, "sample")
        eval_records.append(record)
        if not passed:
            passed_all = False
    return passed_all


def evaluate_length_rules(metrics, args, offsets, length_rules, rule_set_id, eval_records):
    keep_lengths = set()
    for read_len, _offset in offsets:
        passed_all = True
        for rule in length_rules:
            passed, record = evaluate_rule(metrics, args, rule, rule_set_id, "length", read_len)
            eval_records.append(record)
            if not passed:
                passed_all = False
        if passed_all:
            keep_lengths.add(read_len)
    return keep_lengths


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


def write_gate_selection(path, records):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "run_id",
                "sample_id",
                "rule_set_id",
                "tier",
                "selected_for_translon",
                "selected_for_trackhub",
                "selected_reason",
                "pass_lengths_json",
            ]
        )
        writer.writerows(records)


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

    offsets = read_offsets_table(args.best_offset)
    eval_records = []

    apply_for = {item.strip() for item in args.apply_for.split(",") if item.strip()}
    tier_configs = rules.get("track_tiers") or {"good": {}}
    tier_results = {}
    gate_selection_records = []

    for tier_name, tier_config in tier_configs.items():
        tier_rule_set_id = f"{rule_set_id}_{tier_name}"
        sample_rules = resolve_tier_rules(rules, tier_config, "sample")
        length_rules = resolve_tier_rules(rules, tier_config, "length")
        sample_pass = evaluate_sample_rules(metrics, args, sample_rules, tier_rule_set_id, eval_records)
        length_passes = evaluate_length_rules(metrics, args, offsets, length_rules, tier_rule_set_id, eval_records)
        keep_lengths = length_passes if sample_pass else set()
        selected = sample_pass and len(keep_lengths) > 0
        selected_for_translon = "translon" in apply_for and selected
        selected_for_trackhub = "trackhub" in apply_for and selected
        reason = "sample rules pass and at least one passing length" if selected else "rules not met"
        tier_results[tier_name] = {
            "selected_for_translon": selected_for_translon,
            "selected_for_trackhub": selected_for_trackhub,
            "reason": reason,
            "keep_lengths": keep_lengths,
        }
        gate_selection_records.append(
            [
                args.run_id,
                args.sample_id,
                tier_rule_set_id,
                tier_name,
                bool(selected_for_translon),
                bool(selected_for_trackhub),
                reason,
                json.dumps(sorted(keep_lengths)),
            ]
        )

        with open(f"{args.out_prefix}.{tier_name}.pass_lengths.tsv", "w") as handle:
            handle.write("length\n")
            for read_len in sorted(keep_lengths):
                handle.write(f"{read_len}\n")

        write_filtered_offsets(args.best_offset, f"{args.out_prefix}.offsets.{tier_name}.tsv", keep_lengths)

    primary_result = tier_results.get("good", next(iter(tier_results.values())))
    keep_lengths = primary_result["keep_lengths"]
    selected_for_translon = primary_result["selected_for_translon"]
    selected_for_trackhub = primary_result["selected_for_trackhub"]
    reason = primary_result["reason"]

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
    write_gate_selection(f"{args.out_prefix}.gate_selection.tsv", gate_selection_records)
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
                "track_tiers": {
                    tier: {
                        "selected_for_translon": result["selected_for_translon"],
                        "selected_for_trackhub": result["selected_for_trackhub"],
                        "selected_reason": result["reason"],
                        "pass_lengths": sorted(result["keep_lengths"]),
                    }
                    for tier, result in tier_results.items()
                },
            },
            handle,
        )


if __name__ == "__main__":
    main()
