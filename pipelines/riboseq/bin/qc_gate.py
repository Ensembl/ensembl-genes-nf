#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import yaml

from qc_db import connect, upsert_rule_set


OPS = {
    ">=": lambda observed, threshold: observed is not None and observed >= threshold,
    "<=": lambda observed, threshold: observed is not None and observed <= threshold,
    "between": lambda observed, threshold: observed is not None and threshold[0] <= observed <= threshold[1],
}


def get_metric(con, run_id, sample_id, metric, scope, length=None):
    if scope == "sample":
        row = con.execute(
            """
            SELECT value_num
            FROM metrics_core
            WHERE run_id = ? AND sample_id = ? AND metric = ? AND scope = 'sample'
            LIMIT 1
            """,
            [run_id, sample_id, metric],
        ).fetchone()
    else:
        row = con.execute(
            """
            SELECT value_num
            FROM metrics_core
            WHERE run_id = ? AND sample_id = ? AND metric = ? AND scope = 'length' AND length = ?
            LIMIT 1
            """,
            [run_id, sample_id, metric, int(length)],
        ).fetchone()
    return row[0] if row else None


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


def evaluate_rule(con, args, rule, scope, length=None):
    metric = rule["metric"]
    op = rule["op"]
    threshold = rule["threshold"]
    observed = get_metric(con, args.run_id, args.sample_id, metric, scope, length)
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
        None,
        scope,
        int(length) if length is not None else None,
        rule["name"],
        metric,
        op,
        threshold_min,
        threshold_max,
        float(observed) if observed is not None else None,
        None,
        rule["severity"],
        bool(passed),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--best-offset", required=True)
    parser.add_argument("--rules-yaml", required=True)
    parser.add_argument("--rule-set-name", default="default")
    parser.add_argument("--apply-for", default="translon,trackhub")
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    con = connect(args.db)
    rules_text = Path(args.rules_yaml).read_text()
    rules = yaml.safe_load(rules_text) or {}
    rule_set_id = upsert_rule_set(
        con,
        rules.get("rule_set_name", args.rule_set_name),
        int(rules.get("version", 1)),
        rules_text,
    )

    keep_lengths = set()
    eval_records = []

    for read_len, _offset in read_offsets_table(args.best_offset):
        passed_all = True
        for rule in rules.get("length_rules", []):
            passed, record = evaluate_rule(con, args, rule, "length", read_len)
            eval_records.append(record[:2] + (rule_set_id,) + record[3:])
            if not passed:
                passed_all = False
        if passed_all:
            keep_lengths.add(read_len)

    sample_pass = True
    for rule in rules.get("sample_rules", []):
        passed, record = evaluate_rule(con, args, rule, "sample")
        eval_records.append(record[:2] + (rule_set_id,) + record[3:])
        if not passed:
            sample_pass = False

    if eval_records:
        con.executemany(
            """
            INSERT INTO qc_eval(
                eval_id, run_id, sample_id, rule_set_id, scope, length, rule_name, metric, op,
                threshold_num_min, threshold_num_max, observed_num, observed_text, severity, passed, decided_at
            )
            VALUES (uuid(), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, current_timestamp)
            """,
            eval_records,
        )

    apply_for = {item.strip() for item in args.apply_for.split(",") if item.strip()}
    selected = sample_pass and len(keep_lengths) > 0
    selected_for_translon = "translon" in apply_for and selected
    selected_for_trackhub = "trackhub" in apply_for and selected
    reason = "sample rules pass and at least one passing length" if selected else "rules not met"

    con.execute(
        """
        INSERT INTO gate_selection(
            run_id, sample_id, rule_set_id, selected_for_translon,
            selected_for_trackhub, selected_reason, pass_lengths_json
        )
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            args.run_id,
            args.sample_id,
            rule_set_id,
            selected_for_translon,
            selected_for_trackhub,
            reason,
            json.dumps(sorted(keep_lengths)),
        ],
    )

    with open(f"{args.out_prefix}.pass_lengths.tsv", "w") as handle:
        handle.write("length\n")
        for read_len in sorted(keep_lengths):
            handle.write(f"{read_len}\n")

    write_filtered_offsets(args.best_offset, f"{args.out_prefix}.offsets.pass.tsv", keep_lengths)

    with open(f"{args.out_prefix}.qc.json", "w") as handle:
        json.dump(
            {
                "sample_id": args.sample_id,
                "rule_set_id": rule_set_id,
                "selected_for_translon": selected_for_translon,
                "selected_for_trackhub": selected_for_trackhub,
                "pass_lengths": sorted(keep_lengths),
            },
            handle,
        )


if __name__ == "__main__":
    main()
