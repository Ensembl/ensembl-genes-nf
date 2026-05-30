#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path


def fail_or_warn(message: str, action: str, failures: list[str], warnings: list[str]):
    if action in {"fail", "record"}:
        failures.append(message)
    else:
        warnings.append(message)


def main():
    parser = argparse.ArgumentParser(
        description="QC gate for partitioned collapsed FASTA TSV stats"
    )
    parser.add_argument("stats_json", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-total-records", type=int, default=1)
    parser.add_argument("--min-total-counts", type=int, default=1)
    parser.add_argument("--max-catch-all-unique-fraction", type=float, default=0.05)
    parser.add_argument("--max-catch-all-count-fraction", type=float, default=0.05)
    parser.add_argument(
        "--action",
        choices=("fail", "warn", "record"),
        default="fail",
        help=(
            "Whether threshold violations fail the process, only emit warnings, "
            "or record a failed QC report while exiting successfully"
        ),
    )
    args = parser.parse_args()

    stats = json.loads(args.stats_json.read_text())
    sample_id = stats.get("sample_id", args.stats_json.name)
    catch_all = stats.get("catch_all_partition")
    total_records = int(stats.get("total_records", 0))
    total_counts = int(stats.get("total_counts", 0))
    total_unique = int(stats.get("total_unique_sequences", 0))

    failures: list[str] = []
    warnings: list[str] = []

    if total_records < args.min_total_records:
        fail_or_warn(
            f"total_records={total_records} is below {args.min_total_records}",
            args.action,
            failures,
            warnings,
        )
    if total_counts < args.min_total_counts:
        fail_or_warn(
            f"total_counts={total_counts} is below {args.min_total_counts}",
            args.action,
            failures,
            warnings,
        )

    partitions = {
        item.get("partition"): item
        for item in stats.get("partitions", [])
    }
    catch_all_stats = partitions.get(catch_all, {})
    catch_all_unique = int(catch_all_stats.get("unique_sequences", 0))
    catch_all_counts = int(catch_all_stats.get("counts", 0))
    catch_all_unique_fraction = catch_all_unique / total_unique if total_unique else 0.0
    catch_all_count_fraction = catch_all_counts / total_counts if total_counts else 0.0

    if catch_all_unique_fraction > args.max_catch_all_unique_fraction:
        fail_or_warn(
            (
                f"{catch_all} unique fraction {catch_all_unique_fraction:.6f} exceeds "
                f"{args.max_catch_all_unique_fraction:.6f}"
            ),
            args.action,
            failures,
            warnings,
        )
    if catch_all_count_fraction > args.max_catch_all_count_fraction:
        fail_or_warn(
            (
                f"{catch_all} count fraction {catch_all_count_fraction:.6f} exceeds "
                f"{args.max_catch_all_count_fraction:.6f}"
            ),
            args.action,
            failures,
            warnings,
        )

    report = {
        "sample_id": sample_id,
        "passed": not failures,
        "action": args.action,
        "failures": failures,
        "warnings": warnings,
        "observed": {
            "total_records": total_records,
            "total_counts": total_counts,
            "total_unique_sequences": total_unique,
            "catch_all_partition": catch_all,
            "catch_all_unique_sequences": catch_all_unique,
            "catch_all_counts": catch_all_counts,
            "catch_all_unique_fraction": catch_all_unique_fraction,
            "catch_all_count_fraction": catch_all_count_fraction,
        },
        "thresholds": {
            "min_total_records": args.min_total_records,
            "min_total_counts": args.min_total_counts,
            "max_catch_all_unique_fraction": args.max_catch_all_unique_fraction,
            "max_catch_all_count_fraction": args.max_catch_all_count_fraction,
        },
    }

    args.output.write_text(json.dumps(report, indent=2) + "\n")

    for message in warnings:
        print(f"WARNING: {sample_id}: {message}", file=sys.stderr)
    for message in failures:
        print(f"ERROR: {sample_id}: {message}", file=sys.stderr)

    if failures and args.action == "fail":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
