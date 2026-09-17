#!/usr/bin/env python3
"""Validate a stable-ID rewritten GFF3 against its original GFF3 and ID map."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path


ID_ATTRIBUTE_KEYS = {
    "ID",
    "Parent",
    "Name",
    "gene_id",
    "transcript_id",
    "exon_id",
    "protein_id",
}

OWNER_ATTRIBUTE_KEYS = (
    "gene_id",
    "transcript_id",
    "exon_id",
    "protein_id",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a GFF3 produced by rewrite_target_gff3.py."
    )
    parser.add_argument("--original-gff", required=True, type=Path)
    parser.add_argument("--rewritten-gff", required=True, type=Path)
    parser.add_argument("--id-map", required=True, type=Path)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def load_id_map(path: Path) -> dict[str, tuple[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = set(reader.fieldnames or [])

        if {"old_id", "new_id", "new_version"}.issubset(headers):
            old_column = "old_id"
            new_column = "new_id"
            skip_empty_rows = False
        elif {
            "current_stable_id",
            "new_stable_id",
            "new_version",
        }.issubset(headers):
            old_column = "current_stable_id"
            new_column = "new_stable_id"
            skip_empty_rows = True
        else:
            raise ValueError(
                f"{path} must contain either "
                "old_id/new_id/new_version or "
                "current_stable_id/new_stable_id/new_version"
            )

        id_map: dict[str, tuple[str, str]] = {}

        for line_number, row in enumerate(reader, start=2):
            old_id = (row.get(old_column) or "").strip()
            new_id = (row.get(new_column) or "").strip()
            new_version = (row.get("new_version") or "").strip()

            if skip_empty_rows and not old_id and not new_id:
                continue

            if not old_id or not new_id or not new_version:
                raise ValueError(
                    f"{path}:{line_number}: empty old ID, new ID, or version"
                )

            previous = id_map.get(old_id)
            if previous is not None and previous != (new_id, new_version):
                raise ValueError(
                    f"{path}:{line_number}: conflicting mapping for {old_id}"
                )

            id_map[old_id] = (new_id, new_version)

    return id_map


def parse_attributes(text: str) -> list[tuple[str, str]]:
    attributes = []

    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Malformed GFF3 attribute: {item!r}")
        key, value = item.split("=", 1)
        attributes.append((key, value))

    return attributes


def format_attributes(attributes: list[tuple[str, str]]) -> str:
    return ";".join(f"{key}={value}" for key, value in attributes)


def replace_id(value: str, id_map: dict[str, tuple[str, str]]) -> tuple[str, bool]:
    """Replace one ID, preserving an optional GFF3 prefix such as gene:."""
    prefix = ""
    bare_id = value

    if ":" in value:
        prefix, bare_id = value.split(":", 1)
        prefix = f"{prefix}:"

    replacement = id_map.get(bare_id)
    if replacement is None:
        return value, False

    return f"{prefix}{replacement[0]}", replacement[0] != bare_id


def replace_id_list(
    value: str,
    id_map: dict[str, tuple[str, str]],
) -> tuple[str, bool]:
    updated_values = []
    changed = False

    for item in value.split(","):
        updated, item_changed = replace_id(item, id_map)
        updated_values.append(updated)
        changed = changed or item_changed

    return ",".join(updated_values), changed


def expected_attributes(
    original: list[tuple[str, str]],
    id_map: dict[str, tuple[str, str]],
) -> tuple[list[tuple[str, str]], set[str]]:
    owner_id = None
    owner_version = None

    for key, value in original:
        if key not in OWNER_ATTRIBUTE_KEYS:
            continue

        _, bare_id = value.split(":", 1) if ":" in value else ("", value)
        replacement = id_map.get(bare_id)
        if replacement is not None:
            owner_id = bare_id
            owner_version = replacement[1]
            break

    expected = []
    changed_ids: set[str] = set()

    for key, value in original:
        if key in ID_ATTRIBUTE_KEYS:
            updated_value, changed = replace_id_list(value, id_map)
            expected.append((key, updated_value))
            if changed:
                changed_ids.add(key)
        elif key == "version" and owner_id is not None:
            expected.append((key, owner_version))
        else:
            expected.append((key, value))

    return expected, changed_ids


def stable_ids_in_attributes(
    attributes: list[tuple[str, str]],
) -> set[str]:
    stable_ids: set[str] = set()

    for key, value in attributes:
        if key not in ID_ATTRIBUTE_KEYS:
            continue

        for item in value.split(","):
            stable_ids.add(item.split(":", 1)[-1])

    return stable_ids


def main() -> int:
    args = parse_args()
    id_map = load_id_map(args.id_map)

    errors: list[str] = []
    observed_map_ids: set[str] = set()
    changed_features: Counter[str] = Counter()
    original_id_count: Counter[str] = Counter()
    rewritten_id_count: Counter[str] = Counter()
    rewritten_feature_ids: set[str] = set()
    rewritten_parents: list[tuple[int, str]] = []

    with (
        open_text(args.original_gff) as original_handle,
        open_text(args.rewritten_gff) as rewritten_handle,
    ):
        for line_number, (original_line, rewritten_line) in enumerate(
            zip(original_handle, rewritten_handle),
            start=1,
        ):
            if original_line.startswith("#") or rewritten_line.startswith("#"):
                if original_line != rewritten_line:
                    errors.append(
                        f"line {line_number}: comment/header differs"
                    )
                continue

            original_fields = original_line.rstrip("\n").split("\t")
            rewritten_fields = rewritten_line.rstrip("\n").split("\t")

            if len(original_fields) != 9 or len(rewritten_fields) != 9:
                errors.append(
                    f"line {line_number}: expected nine GFF3 columns"
                )
                continue

            if original_fields[:8] != rewritten_fields[:8]:
                errors.append(
                    f"line {line_number}: non-attribute GFF3 columns differ"
                )
                continue

            try:
                original_attributes = parse_attributes(original_fields[8])
                rewritten_attributes = parse_attributes(rewritten_fields[8])
                expected, changed_keys = expected_attributes(
                    original_attributes,
                    id_map,
                )
            except ValueError as error:
                errors.append(f"line {line_number}: {error}")
                continue

            if format_attributes(expected) != format_attributes(rewritten_attributes):
                errors.append(
                    f"line {line_number}: attributes differ from the exact "
                    "changes specified by the ID map"
                )

            feature_type = original_fields[2]
            changed_features[feature_type] += bool(changed_keys)

            for stable_id in stable_ids_in_attributes(original_attributes):
                original_id_count[feature_type] += 1
                if stable_id in id_map:
                    observed_map_ids.add(stable_id)

            for key, value in rewritten_attributes:
                if key == "ID":
                    rewritten_feature_ids.update(value.split(","))

                if key == "Parent":
                    for parent in value.split(","):
                        rewritten_parents.append((line_number, parent))

                if key in ID_ATTRIBUTE_KEYS:
                    rewritten_id_count[feature_type] += len(value.split(","))

        if original_handle.readline() or rewritten_handle.readline():
            errors.append("original and rewritten GFF3 files have different line counts")

    changed_old_ids = {
        old_id
        for old_id, (new_id, _) in id_map.items()
        if old_id != new_id
    }

    old_ids_remaining = set()
    with open_text(args.rewritten_gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            for stable_id in stable_ids_in_attributes(parse_attributes(fields[8])):
                if stable_id in changed_old_ids:
                    old_ids_remaining.add(stable_id)

    if old_ids_remaining:
        errors.append(
            "changed old IDs still present in rewritten GFF3: "
            + ", ".join(sorted(old_ids_remaining)[:10])
        )

    missing_parents = [
        (line_number, parent)
        for line_number, parent in rewritten_parents
        if parent not in rewritten_feature_ids
    ]
    if missing_parents:
        examples = ", ".join(
            f"line {line_number}: {parent}"
            for line_number, parent in missing_parents[:10]
        )
        errors.append(
            f"{len(missing_parents)} rewritten Parent values have no matching "
            f"rewritten ID ({examples})"
        )

    unused_map_ids = set(id_map) - observed_map_ids

    print("GFF3 REWRITE VALIDATION")
    print(f"Original GFF3:  {args.original_gff}")
    print(f"Rewritten GFF3: {args.rewritten_gff}")
    print(f"ID map:         {args.id_map}")
    print(f"Map entries: {len(id_map)}")
    print(f"Map IDs observed in original GFF3: {len(observed_map_ids)}")
    print(f"Map entries not observed in original GFF3: {len(unused_map_ids)}")

    if unused_map_ids:
        print(
            "Examples of unused map IDs: "
            + ", ".join(sorted(unused_map_ids)[:10])
        )

    print("Changed GFF3 features:")
    for feature_type, count in sorted(changed_features.items()):
        if count:
            print(f"  {feature_type}: {count}")

    if errors:
        print(f"\nFAILED: {len(errors)} validation issue(s)", file=sys.stderr)
        for error in errors[:20]:
            print(f"  {error}", file=sys.stderr)
        if len(errors) > 20:
            print(
                f"  ... and {len(errors) - 20} more",
                file=sys.stderr,
            )
        return 1

    print("\nPASSED: rewritten GFF3 differs only by expected stable-ID/version updates.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())