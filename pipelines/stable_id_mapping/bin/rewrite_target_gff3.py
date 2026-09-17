#!/usr/bin/env python3
"""Write a target GFF3 with stable IDs and versions updated from a TSV map."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Write a copy of a target GFF3 with stable IDs and feature "
            "versions updated from a tab-separated ID map."
        )
    )
    parser.add_argument(
        "--input-gff",
        required=True,
        type=Path,
        help="Original target GFF3 file, plain text or gzip-compressed.",
    )
    parser.add_argument(
        "--id-map",
        required=True,
        type=Path,
        help="TSV containing old and new stable IDs plus new_version.",
    )
    parser.add_argument(
        "--output-gff",
        required=True,
        type=Path,
        help="Path for the rewritten plain-text GFF3.",
    )
    return parser.parse_args()


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


def load_id_map(path: Path) -> dict[str, tuple[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])

        if {"old_id", "new_id", "new_version"}.issubset(fieldnames):
            old_id_column = "old_id"
            new_id_column = "new_id"
            skip_empty_rows = False

        elif {
            "current_stable_id",
            "new_stable_id",
            "new_version",
        }.issubset(fieldnames):
            old_id_column = "current_stable_id"
            new_id_column = "new_stable_id"
            skip_empty_rows = True

        else:
            raise ValueError(
                f"{path} must contain either "
                "old_id/new_id/new_version or "
                "current_stable_id/new_stable_id/new_version"
            )

        updates: dict[str, tuple[str, str]] = {}

        for row_number, row in enumerate(reader, start=2):
            old_id = (row.get(old_id_column) or "").strip()
            new_id = (row.get(new_id_column) or "").strip()
            new_version = (row.get("new_version") or "").strip()

            if skip_empty_rows and not old_id and not new_id:
                continue

            if not old_id or not new_id or not new_version:
                raise ValueError(
                    f"{path}:{row_number} has an empty ID or version"
                )

            if old_id in updates:
                raise ValueError(
                    f"{path}:{row_number} repeats old_id {old_id!r}"
                )

            updates[old_id] = (new_id, new_version)

    return updates


def open_gff(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def split_prefix(value: str) -> tuple[str, str]:
    if ":" not in value:
        return "", value

    prefix, stable_id = value.rsplit(":", 1)
    return f"{prefix}:", stable_id


def rewrite_id_value(
    value: str,
    updates: dict[str, tuple[str, str]],
    used_ids: set[str],
) -> tuple[str, int]:
    rewritten_values: list[str] = []
    replacements = 0

    for item in value.split(","):
        prefix, stable_id = split_prefix(item)
        update = updates.get(stable_id)

        if update is None:
            rewritten_values.append(item)
            continue

        new_id, _new_version = update
        rewritten_values.append(f"{prefix}{new_id}")
        used_ids.add(stable_id)
        replacements += 1

    return ",".join(rewritten_values), replacements


def owner_update(
    attributes: list[tuple[str, str | None]],
    updates: dict[str, tuple[str, str]],
) -> tuple[str, str] | None:
    for owner_key in OWNER_ATTRIBUTE_KEYS:
        for key, value in attributes:
            if key != owner_key or value is None:
                continue

            _prefix, stable_id = split_prefix(value)
            update = updates.get(stable_id)

            if update is not None:
                return update

    for key, value in attributes:
        if key != "ID" or value is None:
            continue

        _prefix, stable_id = split_prefix(value)
        update = updates.get(stable_id)

        if update is not None:
            return update

    return None


def rewrite_attributes(
    attributes_text: str,
    updates: dict[str, tuple[str, str]],
    used_ids: set[str],
) -> tuple[str, int, int]:
    attributes: list[tuple[str, str | None]] = []

    for item in attributes_text.split(";"):
        if "=" not in item:
            attributes.append((item, None))
            continue

        key, value = item.split("=", 1)
        attributes.append((key, value))

    update_for_owner = owner_update(attributes, updates)

    rewritten: list[str] = []
    id_replacements = 0
    version_replacements = 0

    for key, value in attributes:
        if value is None:
            rewritten.append(key)
            continue

        if key in ID_ATTRIBUTE_KEYS:
            value, replacements = rewrite_id_value(
                value,
                updates,
                used_ids,
            )
            id_replacements += replacements

        elif key == "version" and update_for_owner is not None:
            _new_id, new_version = update_for_owner
            value = new_version
            version_replacements += 1

        rewritten.append(f"{key}={value}")

    return ";".join(rewritten), id_replacements, version_replacements


def rewrite_gff(
    input_gff: Path,
    output_gff: Path,
    updates: dict[str, tuple[str, str]],
) -> tuple[int, int, int, set[str]]:
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_gff.with_name(f".{output_gff.name}.tmp")

    feature_lines = 0
    id_replacements = 0
    version_replacements = 0
    used_ids: set[str] = set()

    try:
        with (
            open_gff(input_gff) as input_handle,
            temporary_output.open("w", encoding="utf-8") as output_handle,
        ):
            for line in input_handle:
                if line.startswith("#"):
                    output_handle.write(line)
                    continue

                fields = line.rstrip("\n").split("\t")

                if len(fields) != 9:
                    output_handle.write(line)
                    continue

                rewritten_attributes, line_id_replacements, line_version_replacements = (
                    rewrite_attributes(
                        fields[8],
                        updates,
                        used_ids,
                    )
                )

                fields[8] = rewritten_attributes
                output_handle.write("\t".join(fields) + "\n")

                feature_lines += 1
                id_replacements += line_id_replacements
                version_replacements += line_version_replacements

        temporary_output.replace(output_gff)

    except Exception:
        temporary_output.unlink(missing_ok=True)
        raise

    return (
        feature_lines,
        id_replacements,
        version_replacements,
        used_ids,
    )


def main() -> int:
    args = parse_args()
    updates = load_id_map(args.id_map)

    (
        feature_lines,
        id_replacements,
        version_replacements,
        used_ids,
    ) = rewrite_gff(
        args.input_gff,
        args.output_gff,
        updates,
    )

    unused_ids = len(updates) - len(used_ids)

    print(
        f"Feature lines written: {feature_lines}\n"
        f"Stable-ID values rewritten: {id_replacements}\n"
        f"Version attributes rewritten: {version_replacements}\n"
        f"Map IDs not found in GFF3: {unused_ids}",
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
