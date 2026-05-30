#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

from qc_db import connect, insert_artifact, insert_metrics


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def value_or_none(value):
    return None if value == "" else value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--metrics", nargs="*", default=[])
    parser.add_argument("--artifacts", nargs="*", default=[])
    parser.add_argument("--rule-sets", nargs="*", default=[])
    parser.add_argument("--qc-evals", nargs="*", default=[])
    parser.add_argument("--gate-selections", nargs="*", default=[])
    args = parser.parse_args()

    con = connect(args.db)

    metric_rows = []
    for path in args.metrics:
        if not Path(path).exists():
            continue
        for row in read_tsv(path):
            metric_rows.append(
                (
                    row["run_id"],
                    row["sample_id"],
                    row["module"],
                    row["metric"],
                    row["scope"],
                    int(row["length"]) if row["length"] else None,
                    float(row["value_num"]),
                    value_or_none(row["unit"]),
                )
            )
    insert_metrics(con, metric_rows)

    for path in args.artifacts:
        if not Path(path).exists():
            continue
        for row in read_tsv(path):
            insert_artifact(
                con,
                (
                    row["artifact_id"],
                    row["run_id"],
                    row["sample_id"],
                    row["module"],
                    row["name"],
                    row["path"],
                    row["checksum"],
                    int(row["size_bytes"]),
                    row["content_type"],
                ),
            )

    for path in args.rule_sets:
        if not Path(path).exists():
            continue
        for row in read_tsv(path):
            con.execute(
                """
                INSERT INTO qc_rule_set(rule_set_id, name, version, yaml_text, created_at)
                VALUES (?, ?, ?, ?, current_timestamp)
                """,
                [row["rule_set_id"], row["name"], int(row["version"]), row["yaml_text"]],
            )

    for path in args.qc_evals:
        if not Path(path).exists():
            continue
        for row in read_tsv(path):
            con.execute(
                """
                INSERT INTO qc_eval(
                    eval_id, run_id, sample_id, rule_set_id, scope, length, rule_name, metric, op,
                    threshold_num_min, threshold_num_max, observed_num, observed_text, severity, passed, decided_at
                )
                VALUES (uuid(), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, current_timestamp)
                """,
                [
                    row["run_id"],
                    row["sample_id"],
                    row["rule_set_id"],
                    row["scope"],
                    int(row["length"]) if row["length"] else None,
                    row["rule_name"],
                    row["metric"],
                    row["op"],
                    float(row["threshold_num_min"]),
                    float(row["threshold_num_max"]) if row["threshold_num_max"] else None,
                    float(row["observed_num"]) if row["observed_num"] else None,
                    value_or_none(row["observed_text"]),
                    row["severity"],
                    row["passed"].lower() == "true",
                ],
            )

    for path in args.gate_selections:
        if not Path(path).exists():
            continue
        for row in read_tsv(path):
            con.execute(
                """
                INSERT INTO gate_selection(
                    run_id, sample_id, rule_set_id, selected_for_translon,
                    selected_for_trackhub, selected_reason, pass_lengths_json
                )
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    row["run_id"],
                    row["sample_id"],
                    row["rule_set_id"],
                    row["selected_for_translon"].lower() == "true",
                    row["selected_for_trackhub"].lower() == "true",
                    row["selected_reason"],
                    row["pass_lengths_json"],
                ],
            )


if __name__ == "__main__":
    main()
