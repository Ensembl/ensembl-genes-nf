#!/usr/bin/env python3
"""Validate and normalize the small, explicit long-read input manifest."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

REQUIRED = {"id", "run_id", "cohort_id", "technology", "assembly_id", "read_set_id", "sample_id", "library_id", "read_representation", "source_uri"}
BAM_REPRESENTATIONS = {"bam", "ccs", "raw_subreads"}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonicalize_approved_row(row: dict[str, str], cohort_id: str) -> dict[str, str]:
    """Map the inventory product contract onto the processing manifest contract."""
    row = dict(row)
    run_id = row.get("id") or row.get("run_accession") or row.get("accession") or row.get("sample_id")
    if not run_id:
        raise SystemExit("approved manifest row has no run_accession, accession, id, or sample_id")
    classification = (row.get("classification") or "").strip().upper()
    representation = (row.get("read_representation") or "").strip().lower()
    is_bam = bool(row.get("bam") or row.get("ccs", "").lower().endswith(".bam") or classification.endswith("_BAM"))
    if is_bam:
        if representation in {"raw_subread", "raw_subreads"} or "SUBREAD" in classification:
            representation = "raw_subreads"
        elif representation == "ccs" or "CCS" in classification:
            representation = "ccs"
        else:
            representation = "bam"
    else:
        representation = "fastq"
    source_uri = row.get("source_uri") or row.get("selected_artifact_uri") or row.get("input", "")
    row.update(
        {
            "id": run_id,
            "run_id": row.get("run_id") or run_id,
            "cohort_id": row.get("cohort_id") or cohort_id,
            "technology": row.get("technology") or "unknown",
            "assembly_id": row.get("assembly_id") or "reference",
            "read_set_id": row.get("read_set_id") or run_id,
            "sample_id": row.get("sample_id") or run_id,
            "library_id": row.get("library_id") or run_id,
            "read_representation": representation,
            "source_uri": source_uri,
            "source_read_representation": row.get("read_representation", ""),
        }
    )
    if not row.get("input"):
        row["input"] = source_uri
    return row


def cached_input(row: dict[str, str], cache_dir: Path | None) -> Path | None:
    if cache_dir is None or not row.get("selected_artifact_md5"):
        return None
    classification = (row.get("classification") or "").strip()
    basename = row.get("selected_artifact_basename") or Path(row.get("input", "")).name
    if not classification or not basename:
        return None
    candidate = cache_dir / classification / row["selected_artifact_md5"] / basename
    return candidate if candidate.is_file() else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=["fastq", "bam", "approved_manifest", "retrospective"])
    parser.add_argument("--manifest-root", type=Path)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--ccs-mode", choices=["disabled", "pbccs"], default="disabled")
    parser.add_argument("--allow-unresolved-inputs", action="store_true")
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--normalized", required=True, type=Path)
    parser.add_argument("--checksums", required=True, type=Path)
    parser.add_argument("--manifest-checksum", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    args = parser.parse_args()
    if not args.manifest.is_file():
        raise SystemExit(f"manifest does not exist: {args.manifest}")
    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit("manifest contains no data rows")
    if args.mode == "approved_manifest":
        rows = [canonicalize_approved_row(row, args.cohort_id) for row in rows]
    missing = sorted(REQUIRED - set(rows[0]))
    if missing:
        raise SystemExit("manifest is missing required columns: " + ", ".join(missing))
    args.normalized.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    if "input" not in fieldnames:
        candidates = [name for name in ("fastq", "bam", "ccs", "path") if name in fieldnames]
        if len(candidates) != 1 and not (args.allow_unresolved_inputs and "source_uri" in fieldnames):
            raise SystemExit("manifest must contain exactly one input, fastq, bam, ccs, or path column")
        fieldnames.append("input")
    resolved = []
    for row in rows:
        row = dict(row)
        if not row.get("cohort_id"):
            row["cohort_id"] = args.cohort_id
        if row["cohort_id"] != args.cohort_id:
            raise SystemExit(f"row {row.get('id')} belongs to cohort {row['cohort_id']!r}, expected {args.cohort_id!r}")
        cached_path = cached_input(row, args.cache_dir) if args.mode == "approved_manifest" else None
        if cached_path:
            row["input"] = str(cached_path)
        if not row.get("input"):
            candidates = [row.get(name) for name in ("fastq", "bam", "ccs", "path") if row.get(name)]
            if not candidates and args.allow_unresolved_inputs and row.get("source_uri"):
                candidates = [row["source_uri"]]
            if len(candidates) != 1:
                raise SystemExit(f"row {row.get('id')} must identify one input path")
            row["input"] = candidates[0]
        input_value = row["input"]
        input_path = Path(input_value)
        if not input_path.is_absolute():
            root = args.manifest_root or args.manifest.parent
            input_path = (root / input_path).resolve()
        if not input_path.is_file():
            source_uri = row.get("source_uri", "")
            remote_source = any(source_uri.startswith(prefix) for prefix in ("http://", "https://", "ftp://", "s3://", "gs://"))
            if not (args.allow_unresolved_inputs and remote_source):
                raise SystemExit(f"input for {row.get('id')} does not exist: {input_path}")
            row["input"] = source_uri
            row["input_sha256"] = ""
        else:
            row["input"] = str(input_path)
            row["input_sha256"] = sha256(input_path)
            expected_md5 = row.get("selected_artifact_md5", "").strip().lower()
            if expected_md5 and len(expected_md5) == 32 and md5(input_path) != expected_md5:
                raise SystemExit(f"input checksum mismatch for {row.get('id')}: expected MD5 {expected_md5}")
        if args.mode == "retrospective":
            model_columns = ("bed12", "model_bed12", "model", "bed12_path")
            for model_column in model_columns:
                if not row.get(model_column):
                    continue
                model_path = Path(row[model_column])
                if not model_path.is_absolute():
                    root = args.manifest_root or args.manifest.parent
                    model_path = (root / model_path).resolve()
                if not model_path.is_file():
                    raise SystemExit(f"model artifact for {row.get('id')} does not exist: {model_path}")
                row[model_column] = str(model_path)
                row[f"{model_column}_sha256"] = sha256(model_path)
        representation = row.get("read_representation")
        if representation not in {"fastq", *BAM_REPRESENTATIONS}:
            raise SystemExit(f"row {row.get('id')} has unsupported read_representation: {representation}")
        if representation == "raw_subreads" and args.ccs_mode != "pbccs":
            raise SystemExit(f"row {row.get('id')} contains raw_subreads but --ccs-mode is not pbccs")
        if args.mode in {"bam", "approved_manifest", "retrospective"} and representation in BAM_REPRESENTATIONS and input_path.is_file():
            suffixes = [".pbi", ".csi", ".bai"] if representation == "raw_subreads" and args.ccs_mode == "pbccs" else [".csi", ".bai"]
            candidates = [input_path.with_name(input_path.name + suffix) for suffix in suffixes]
            candidates.extend(input_path.with_suffix(suffix) for suffix in suffixes)
            index_path = next((candidate for candidate in candidates if candidate.is_file()), None)
            if index_path is None:
                expected = " or ".join(str(candidate) for candidate in candidates)
                raise SystemExit(f"indexed BAM sidecar is required for {row.get('id')}: expected {expected}")
            row["input_index"] = str(index_path.resolve())
            row["input_index_sha256"] = sha256(index_path)
        row["manifest_mode"] = args.mode
        resolved.append(row)
        fieldnames.extend(name for name in ("input_sha256", "input_index", "input_index_sha256", "manifest_mode") if name not in fieldnames)
        if args.mode == "retrospective":
            for model_column in ("bed12", "model_bed12", "model", "bed12_path"):
                checksum_column = f"{model_column}_sha256"
                if row.get(model_column) and checksum_column not in fieldnames:
                    fieldnames.append(checksum_column)
    fieldnames = list(dict.fromkeys(fieldnames))
    with args.normalized.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(resolved)
    with args.checksums.open("w") as handle:
        handle.write("id\tpath\tsha256\n")
        for row in resolved:
            handle.write(f"{row['id']}\t{row['input']}\t{row['input_sha256']}\n")
            if row.get("input_index"):
                handle.write(f"{row['id']}[index]\t{row['input_index']}\t{row['input_index_sha256']}\n")
    manifest_digest = sha256(args.normalized)
    args.manifest_checksum.write_text(f"{manifest_digest}  {args.normalized.name}\n")
    with args.status.open("w") as handle:
        handle.write("id\tstatus\tmessage\n")
        for row in resolved:
            handle.write(f"{row['id']}\tVALID\tinput checksum recorded\n")
    args.status.with_suffix(".json").write_text(json.dumps({"rows": len(resolved), "mode": args.mode}, indent=2) + "\n")


if __name__ == "__main__":
    main()
