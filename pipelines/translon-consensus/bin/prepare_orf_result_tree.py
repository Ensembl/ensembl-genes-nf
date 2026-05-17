#!/usr/bin/env python3
"""Stage messy translon caller outputs into a consistent tool/sample directory tree."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
from pathlib import Path

from translon_db_standardise import detect_parser, normalise_sample


SUPPORTED_SUFFIXES = {".bed", ".bed12", ".gtf", ".gff", ".gff3", ".csv", ".tsv", ".txt"}


def clean_token(value: str) -> str:
    value = re.sub(r"\.(known|novel)$", "", value)
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return value.strip("_") or "unknown"


def compressed_suffix(path: Path) -> str:
    suffixes = path.suffixes
    if len(suffixes) >= 2 and suffixes[-1] == ".gz":
        return "".join(suffixes[-2:])
    return path.suffix


def infer_sample(path: Path) -> str:
    name = path.name
    for pattern in [
        r"^ncORFs_(.+)\.bed(?:\.gz)?$",
        r"^(PRICE_Annotations_|RiboTIE_Annotations_|ORFQuant_Annotations_|iRibo_Annotations_)?(.+)\.(bigBed|bed12|bed|gtf|gff3?|csv|tsv)(?:\.gz)?$",
        r"^(.+)\.(known|novel)\.(bed12|bed)(?:\.gz)?$",
    ]:
        match = re.match(pattern, name)
        if match:
            candidate = next(group for group in reversed(match.groups()) if group and group not in {"known", "novel", "bed", "bed12", "bigBed", "gtf", "gff", "gff3", "csv", "tsv"})
            return normalise_sample(candidate)
    return normalise_sample(path.stem)


def iter_file_list(path: Path) -> list[dict[str, str]]:
    rows = []
    with path.open() as handle:
        first = handle.readline()
        handle.seek(0)
        if "\t" in first and "path" in first:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                rows.append({key: row.get(key, "") for key in ["path", "tool", "sample_id"]})
        else:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    rows.append({"path": line, "tool": "", "sample_id": ""})
    return rows


def discover(input_root: Path | None, file_list: Path | None) -> list[dict[str, str]]:
    if file_list:
        return iter_file_list(file_list)
    assert input_root is not None
    rows = []
    for path in sorted(input_root.rglob("*")):
        if not path.is_file():
            continue
        suffix = compressed_suffix(path).replace(".gz", "")
        if suffix in SUPPORTED_SUFFIXES:
            rows.append({"path": str(path), "tool": "", "sample_id": ""})
    return rows


def stage_file(src: Path, dest: Path, copy: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if copy:
        shutil.copy2(src, dest)
    else:
        os.symlink(src.resolve(), dest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--file-list", type=Path, help="Plain path list or TSV with path/tool/sample_id columns")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--copy", action="store_true", help="Copy files instead of symlinking")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.input_root and not args.file_list:
        raise SystemExit("Provide --input-root or --file-list")

    clean_root = args.out_dir / "clean_orf_inputs"
    clean_root.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    seen: dict[Path, int] = {}

    for row in discover(args.input_root, args.file_list):
        src = Path(row["path"])
        if not src.exists():
            manifest_rows.append({"original_path": str(src), "stage_status": "missing"})
            continue
        tool_hint = row.get("tool") or src.parent.name
        status, parser_name, detected_tool = detect_parser(src, tool_hint)
        tool = clean_token(row.get("tool") or detected_tool or tool_hint)
        sample_id = clean_token(row.get("sample_id") or infer_sample(src))
        suffix = compressed_suffix(src)
        dest = clean_root / tool / f"{sample_id}{suffix}"
        if dest in seen:
            seen[dest] += 1
            dest = clean_root / tool / f"{sample_id}.{seen[dest]}{suffix}"
        else:
            seen[dest] = 1
        if status == "ok":
            stage_file(src, dest, args.copy)
            stage_status = "staged_copy" if args.copy else "staged_symlink"
        else:
            stage_status = "unsupported"
        manifest_rows.append(
            {
                "original_path": str(src),
                "clean_path": str(dest) if status == "ok" else "",
                "tool": tool,
                "sample_id": sample_id,
                "parser_status": status,
                "parser_name": parser_name,
                "detected_tool": detected_tool,
                "stage_status": stage_status,
            }
        )

    manifest = args.out_dir / "clean_orf_inputs_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        fieldnames = [
            "original_path",
            "clean_path",
            "tool",
            "sample_id",
            "parser_status",
            "parser_name",
            "detected_tool",
            "stage_status",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {clean_root}")
    print(f"Wrote {manifest}")


if __name__ == "__main__":
    main()
