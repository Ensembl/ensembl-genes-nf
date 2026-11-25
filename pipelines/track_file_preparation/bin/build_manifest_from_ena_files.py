#!/usr/bin/env python3
"""Build a track manifest from an ENA submission files.tsv inventory.

ENA files.tsv is treated as the authoritative BAM/run list. Companion GTF and
STAR junction paths are resolved from their explicitly supplied roots.
Diagnostics are always written; in strict mode the manifest is withheld if
any requested input is missing or ambiguous.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


GCA_RE = re.compile(r"GCA_[0-9]+\.[0-9]+")
MANIFEST_FIELDS = [
    "entity_id", "entity_type", "gca_accession", "assembly_release", "bam", "gtf",
    "track_types", "bai", "sj_out_tab", "source_entity_ids", "run_accession", "sample_id",
]
AUDIT_FIELDS = ["entity_id", "gca_accession", "bam", "gtf", "sj_out_tab", "status", "detail"]
EXCEPTION_FIELDS = ["entity_id", "exception_class", "path", "detail"]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files-tsv", required=True, type=Path,
                        help="ENA submission files.tsv; file_path/run_accession are required")
    parser.add_argument("--gca", required=True, help="Assembly accession, e.g. GCA_052040795.2")
    parser.add_argument("--assembly-release", required=True)
    parser.add_argument("--gtf-dir", required=True, type=Path,
                        help="Directory containing <run>_Aligned.sortedByCoord.out.gtf")
    parser.add_argument("--sj-dir", required=True, type=Path,
                        help="STAR output root containing <run>__STARpass1/SJ.out.tab")
    parser.add_argument("--track-types", default="coverage,gene_model,splice_junction")
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--no-strict", dest="strict", action="store_false",
                        help="Write available rows despite missing requested companions")
    parser.set_defaults(strict=True)
    return parser.parse_args(argv)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def one_file(path: Path, pattern: str) -> Tuple[Optional[Path], str]:
    matches = sorted(path.glob(pattern))
    if len(matches) == 1:
        return matches[0].resolve(), ""
    if not matches:
        return None, f"no match for {pattern}"
    return None, f"multiple matches for {pattern}: {', '.join(map(str, matches))}"


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if not GCA_RE.fullmatch(args.gca):
        print(f"invalid --gca: {args.gca}", file=sys.stderr)
        return 2
    tracks = [part.strip() for part in args.track_types.split(",") if part.strip()]
    allowed = {"coverage", "gene_model", "splice_junction"}
    if not tracks or set(tracks) - allowed:
        print("--track-types must be a comma-separated subset of: " + ", ".join(sorted(allowed)),
              file=sys.stderr)
        return 2
    tracks = [track for track in ("coverage", "gene_model", "splice_junction") if track in tracks]

    try:
        with args.files_tsv.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"file_path", "file_type", "run_accession"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError("files.tsv missing columns: " + ", ".join(sorted(missing)))
            input_rows = [row for row in reader if row.get("file_type", "").lower() == "bam"]
    except (OSError, ValueError, csv.Error) as exc:
        print(str(exc), file=sys.stderr)
        return 2

    rows: List[Dict[str, str]] = []
    audit: List[Dict[str, str]] = []
    exceptions: List[Dict[str, str]] = []
    seen = set()
    for source in input_rows:
        run = (source.get("run_accession") or "").strip()
        bam_text = (source.get("file_path") or "").strip()
        sample = (source.get("sample_accession") or "").strip()
        problems: List[str] = []
        if not run:
            exceptions.append({"entity_id": "", "exception_class": "missing_run_accession",
                               "path": bam_text, "detail": "BAM row has empty run_accession"})
            continue
        if run in seen:
            exceptions.append({"entity_id": run, "exception_class": "duplicate_run",
                               "path": bam_text, "detail": "run occurs more than once in files.tsv"})
            continue
        seen.add(run)
        bam = Path(bam_text)
        if not bam.is_file():
            problems.append(f"BAM does not exist: {bam}")
        gtf, gtf_problem = one_file(args.gtf_dir, f"{run}_Aligned.sortedByCoord.out.gtf")
        sj, sj_problem = one_file(args.sj_dir, f"{run}__STARpass1/SJ.out.tab")
        if "gene_model" in tracks and gtf is None:
            problems.append("GTF " + gtf_problem)
        if "splice_junction" in tracks and sj is None:
            problems.append("STAR junction " + sj_problem)
        bai = Path(str(bam) + ".bai")
        csi = Path(str(bam) + ".csi")
        index = bai if bai.is_file() else (csi if csi.is_file() else None)
        rows.append({
            "entity_id": run, "entity_type": "run", "gca_accession": args.gca,
            "assembly_release": args.assembly_release, "bam": str(bam.resolve()) if bam.is_file() else bam_text,
            "gtf": str(gtf) if gtf else "", "track_types": ",".join(tracks),
            "bai": str(index.resolve()) if index else "", "sj_out_tab": str(sj) if sj else "",
            "source_entity_ids": "", "run_accession": run, "sample_id": sample,
        })
        audit.append({"entity_id": run, "gca_accession": args.gca, "bam": bam_text,
                      "gtf": str(gtf) if gtf else "", "sj_out_tab": str(sj) if sj else "",
                      "status": "resolved" if not problems else "incomplete", "detail": "; ".join(problems)})
        for problem in problems:
            kind = "missing_bam" if problem.startswith("BAM") else (
                "missing_gtf" if problem.startswith("GTF") else "missing_sj")
            exceptions.append({"entity_id": run, "exception_class": kind,
                               "path": bam_text, "detail": problem})

    output = args.output_manifest
    write_tsv(output.with_name("manifest_audit.tsv"), AUDIT_FIELDS, audit)
    write_tsv(output.with_name("manifest_exceptions.tsv"), EXCEPTION_FIELDS, exceptions)
    can_publish = bool(rows) and (not args.strict or not exceptions)
    output.with_name("manifest_summary.json").write_text(json.dumps({
        "source_files_tsv": str(args.files_tsv.resolve()), "gca_accession": args.gca,
        "assembly_release": args.assembly_release, "input_bam_rows": len(input_rows),
        "manifest_rows": len(rows), "audit_rows": len(audit), "exception_rows": len(exceptions),
        "track_types": tracks, "strict": args.strict, "published": can_publish,
    }, indent=2, sort_keys=True) + "\n")
    if can_publish:
        write_tsv(output, MANIFEST_FIELDS, rows)
    if exceptions and args.strict:
        return 2
    if not rows:
        print("No BAM rows found in files.tsv", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
