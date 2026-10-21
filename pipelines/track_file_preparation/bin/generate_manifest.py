#!/usr/bin/env python3
"""Discover run-level track inputs for track_file_preparation.

The script deliberately emits only ``run`` entities.  Sample and merged
entities require an explicit grouping step and are outside discovery.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


MANIFEST_FIELDS = [
    "entity_id", "entity_type", "gca_accession", "assembly_release", "bam", "gtf",
    "track_types", "bai", "sj_out_tab", "source_entity_ids", "run_accession", "sample_id",
]
AUDIT_FIELDS = [
    "source_id", "source_type", "entity_id", "entity_type", "gca_accession", "source_path",
    "bam", "bai", "gtf", "sj_out_tab", "matching_rule", "status", "warning",
]
EXCEPTION_FIELDS = ["source_id", "source_type", "entity_id", "exception_class", "path", "detail"]
GCA_RE = re.compile(r"GCA_[0-9]+\.[0-9]+")
RUN_RE = re.compile(r"(?:ERR|SRR|DRR)[0-9]+")
BAM_SUFFIX = "_Aligned.sortedByCoord.out.bam"
GTF_SUFFIXES = (
    "_Aligned.sortedByCoord.out.gtf",
    ".scallop.gtf",
    ".stringtie.gtf",
    ".gtf",
)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=("anno", "gdm", "main"))
    parser.add_argument("--base-dir", required=True, type=Path)
    parser.add_argument("--entity", help="Optional run ID or directory below base-dir")
    parser.add_argument("--annotation-dir", type=Path, help="Override the mode's annotation directory")
    parser.add_argument("--gca", help="GCA accession; inferred from the base path when omitted")
    parser.add_argument("--assembly-release", help="Assembly release; defaults to --gca")
    parser.add_argument("--gdm-duplicate-policy", choices=("error", "newest"), default="error",
                        help="How to resolve repeated GDM run BAMs; newest selects the unique latest mtime")
    parser.add_argument(
        "--track-types",
        help="Comma-delimited tracks; defaults to all for main/anno and coverage,splice_junction for gdm",
    )
    parser.add_argument("--output-manifest", type=Path, required=True)
    parser.add_argument("--output-audit", type=Path)
    parser.add_argument("--output-exceptions", type=Path)
    parser.add_argument("--output-summary", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-strict", dest="strict", action="store_false")
    parser.set_defaults(strict=True)
    return parser.parse_args(argv)


def normalise_tracks(raw: Optional[str], mode: str) -> List[str]:
    if raw is None or not raw.strip():
        raw = "coverage,splice_junction" if mode == "gdm" else "all"
    selected = [part.strip() for part in raw.split(",") if part.strip()]
    if not selected:
        selected = (
            ["coverage", "splice_junction"]
            if mode == "gdm"
            else ["coverage", "gene_model", "splice_junction"]
        )
    if selected == ["all"]:
        return ["coverage", "gene_model", "splice_junction"]
    allowed = {"coverage", "gene_model", "splice_junction"}
    invalid = sorted(set(selected) - allowed)
    if invalid:
        raise ValueError("unsupported track_types: " + ", ".join(invalid))
    return [name for name in ("coverage", "gene_model", "splice_junction") if name in selected]


def absolute(path: Path) -> Path:
    return Path(os.path.abspath(os.fspath(path)))


def readable_dir(path: Path) -> bool:
    return path.is_dir() and os.access(path, os.R_OK | os.X_OK)


def infer_gca(path: Path) -> Optional[str]:
    path_text = str(path)
    match = GCA_RE.search(path_text)
    if match:
        return match.group(0)
    # GDM species directories encode the assembly as e.g. gca000003025v6.
    slug_match = re.search(r"(?i)gca[_-]?([0-9]{9,})v([0-9]+)", path_text)
    if slug_match:
        return f"GCA_{slug_match.group(1)}.{slug_match.group(2)}"
    return None


def run_id_from_bam(path: Path) -> str:
    if not path.name.endswith(BAM_SUFFIX):
        raise ValueError("unsupported BAM filename")
    return path.name[: -len(BAM_SUFFIX)]


def find_gtf(annotation_dir: Optional[Path], run_id: str) -> Optional[Path]:
    if annotation_dir is None:
        return None
    matches = [annotation_dir / (run_id + suffix) for suffix in GTF_SUFFIXES]
    existing = [path for path in matches if path.is_file()]
    if len(existing) > 1:
        raise ValueError("ambiguous GTF match: " + ", ".join(map(str, existing)))
    return existing[0] if existing else None


def companion_paths(bam: Path, annotation_dir: Optional[Path], run_id: str) -> Dict[str, Optional[Path]]:
    bai = Path(str(bam) + ".bai")
    csi = Path(str(bam) + ".csi")
    # Prefer the broadly supported BAI when both valid index formats are present.
    index = bai if bai.is_file() else (csi if csi.is_file() else None)
    sj = bam.parent / (run_id + "_SJ.out.tab")
    return {"bai": index, "sj_out_tab": sj if sj.is_file() else None,
            "gtf": find_gtf(annotation_dir, run_id)}


def mode_roots(args: argparse.Namespace) -> Tuple[Path, Optional[Path]]:
    base = absolute(args.base_dir)
    if args.entity and (base / args.entity).is_dir():
        base = base / args.entity
    if args.annotation_dir:
        annotation = absolute(args.annotation_dir)
    elif args.mode == "anno":
        annotation = base / "scallop_output"
    elif args.mode == "main":
        annotation = next((base / "rnaseq" / name / "scallop" for name in ("output", "outputs")
                           if (base / "rnaseq" / name / "scallop").is_dir()), None)
    else:
        annotation = None
    return base, annotation


def main_alignment_root(base: Path) -> Path:
    candidates = [base / "rnaseq" / name for name in ("output", "outputs")]
    for candidate in candidates:
        if any(candidate.glob("*" + BAM_SUFFIX)):
            return candidate
    return next((candidate for candidate in candidates if candidate.is_dir()), base)


def candidate_bams(args: argparse.Namespace, base: Path) -> Iterable[Path]:
    if args.mode == "anno":
        root = base / "star_output" if (base / "star_output").is_dir() else base
        yield from sorted(root.glob("*" + BAM_SUFFIX))
    elif args.mode == "main":
        root = main_alignment_root(base)
        yield from sorted(root.glob("*" + BAM_SUFFIX))
    else:
        for path in sorted(base.rglob("*" + BAM_SUFFIX)):
            if path.parent.name == "alignment":
                yield path


def entity_matches(args: argparse.Namespace, run_id: str, path: Path) -> bool:
    if not args.entity or (absolute(args.base_dir) / args.entity).is_dir():
        return True
    return run_id == args.entity


def validate_gdm_path(path: Path) -> Optional[str]:
    alignment = path.parent
    run_dir = alignment.parent
    if alignment.name != "alignment" or not run_dir.name:
        return "expected <assembly>/<optional-platform>/<run>/alignment"
    try:
        run_id = run_id_from_bam(path)
    except ValueError:
        return "invalid BAM name"
    if run_dir.name != run_id:
        return "run directory and BAM basename disagree"
    if not RUN_RE.fullmatch(run_id):
        return "run directory is not an ERR/SRR/DRR accession"
    return None


def collapse_gdm_candidates(args: argparse.Namespace, candidates: Sequence[Path]) -> Tuple[List[Path], Dict[Path, str]]:
    """Optionally select one uniquely newest path for duplicate GDM run IDs."""
    if args.mode != "gdm" or args.gdm_duplicate_policy != "newest":
        return list(candidates), {}

    groups: Dict[Tuple[str, str], List[Path]] = {}
    passthrough: List[Path] = []
    for bam in candidates:
        run_id = run_id_from_bam(bam)
        gca = args.gca or infer_gca(bam)
        if not gca or not GCA_RE.fullmatch(gca) or validate_gdm_path(bam):
            passthrough.append(bam)
            continue
        groups.setdefault((gca, run_id), []).append(bam)

    selected = list(passthrough)
    warnings: Dict[Path, str] = {}
    for key, paths in groups.items():
        if len(paths) == 1:
            selected.extend(paths)
            continue
        ordered = sorted(paths, key=lambda path: (path.stat().st_mtime_ns, str(path)))
        newest = ordered[-1]
        if len(ordered) > 1 and ordered[-2].stat().st_mtime_ns == newest.stat().st_mtime_ns:
            # A timestamp tie stays unresolved and will be reported by duplicate checks.
            selected.extend(paths)
            continue
        selected.append(newest)
        warnings[newest] = f"selected newest of {len(paths)} candidate BAM paths by mtime for {key[1]}"
    return sorted(selected), warnings


def make_row(args: argparse.Namespace, bam: Optional[Path], gtf: Optional[Path], sj: Optional[Path],
             bai: Optional[Path], run_id: str, gca: str) -> Dict[str, str]:
    tracks = normalise_tracks(args.track_types, args.mode)
    required = {"coverage": bam, "gene_model": gtf, "splice_junction": sj}
    missing = [name for name in tracks if required[name] is None]
    if missing:
        raise ValueError("missing required " + ", ".join(missing))
    values = {
        "entity_id": run_id, "entity_type": "run", "gca_accession": gca,
        "assembly_release": args.assembly_release or gca,
        "bam": str(absolute(bam)) if bam else "", "gtf": str(absolute(gtf)) if gtf else "",
        "track_types": ",".join(tracks), "bai": str(absolute(bai)) if bai else "",
        "sj_out_tab": str(absolute(sj)) if sj else "", "source_entity_ids": "",
        "run_accession": run_id if RUN_RE.fullmatch(run_id) else "", "sample_id": "",
    }
    for key, value in values.items():
        if "\t" in value or "\n" in value or (key in ("bam", "gtf", "bai", "sj_out_tab") and "," in value):
            raise ValueError("invalid delimiter in " + key)
    return values


def discover(args: argparse.Namespace) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], List[Dict[str, str]]]:
    base, annotation = mode_roots(args)
    exceptions: List[Dict[str, str]] = []
    audit: List[Dict[str, str]] = []
    rows: List[Dict[str, str]] = []
    if not readable_dir(base):
        exceptions.append({"source_id": args.entity or base.name, "source_type": args.mode,
                           "entity_id": "", "exception_class": "unreadable_root", "path": str(base),
                           "detail": "directory does not exist or is not readable"})
        return rows, audit, exceptions
    base_gca = args.gca or (None if args.mode == "gdm" else infer_gca(base))
    if args.mode != "gdm" and (not base_gca or not GCA_RE.fullmatch(base_gca)):
        exceptions.append({"source_id": args.entity or base.name, "source_type": args.mode,
                           "entity_id": "", "exception_class": "assembly_mismatch", "path": str(base),
                           "detail": "provide --gca or include a GCA accession in the base path"})
        return rows, audit, exceptions
    seen: Dict[Tuple[str, str], Path] = {}
    candidates, duplicate_warnings = collapse_gdm_candidates(args, list(candidate_bams(args, base)))
    for bam in candidates:
        try:
            run_id = run_id_from_bam(bam)
            if not entity_matches(args, run_id, bam):
                continue
            path_gca = infer_gca(bam)
            gca = (args.gca or path_gca) if args.mode == "gdm" else base_gca
            if not gca or not GCA_RE.fullmatch(gca):
                raise ValueError("could not infer a valid GCA accession from input path; provide --gca")
            if path_gca and path_gca != gca:
                raise ValueError("GCA in input path does not match --gca or selected assembly")
            if args.mode == "gdm":
                problem = validate_gdm_path(bam)
                if problem:
                    raise ValueError(problem)
            key = (gca, run_id)
            if key in seen:
                raise ValueError("duplicate or ambiguous run: " + str(seen[key]))
            seen[key] = bam
            companions = companion_paths(bam, annotation, run_id)
            row = make_row(args, bam, companions["gtf"], companions["sj_out_tab"], companions["bai"], run_id, gca)
            rows.append(row)
            audit.append({"source_id": args.entity or (gca if args.mode == "gdm" else base.name), "source_type": args.mode,
                          "entity_id": run_id, "entity_type": "run", "gca_accession": gca,
                          "source_path": str(absolute(bam.parent)), "bam": row["bam"], "bai": row["bai"],
                          "gtf": row["gtf"], "sj_out_tab": row["sj_out_tab"],
                          "matching_rule": args.mode, "status": "resolved",
                          "warning": duplicate_warnings.get(bam, "")})
        except (OSError, ValueError) as exc:
            detail = str(exc)
            if detail.startswith("missing required "):
                missing_name = detail[len("missing required "):].split(",", 1)[0]
                exception_class = "missing_" + missing_name
            elif "duplicate or ambiguous" in detail:
                exception_class = "duplicate_entity"
            elif "ambiguous" in detail:
                exception_class = "ambiguous_match"
            elif "GCA" in detail:
                exception_class = "assembly_mismatch"
            else:
                exception_class = "unsupported_layout"
            exceptions.append({"source_id": args.entity or base.name, "source_type": args.mode,
                               "entity_id": locals().get("run_id", ""), "exception_class": exception_class,
                               "path": str(bam), "detail": detail})
    if not rows and not exceptions:
        exceptions.append({"source_id": args.entity or base.name, "source_type": args.mode,
                           "entity_id": "", "exception_class": "unsupported_layout", "path": str(base),
                           "detail": "no matching STAR BAM files found"})
    rows.sort(key=lambda row: (row["gca_accession"], row["entity_id"]))
    audit.sort(key=lambda row: (row["gca_accession"], row["entity_id"]))
    return rows, audit, exceptions


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        normalise_tracks(args.track_types, args.mode)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    rows, audit, exceptions = discover(args)
    output_audit = args.output_audit or args.output_manifest.with_name("manifest_audit.tsv")
    output_exceptions = args.output_exceptions or args.output_manifest.with_name("manifest_exceptions.tsv")
    output_summary = args.output_summary or args.output_manifest.with_name("manifest_summary.json")
    write_tsv(output_audit, AUDIT_FIELDS, audit)
    write_tsv(output_exceptions, EXCEPTION_FIELDS, exceptions)
    summary = {"mode": args.mode, "base_dir": str(absolute(args.base_dir)),
               "gca_accessions": sorted({row["gca_accession"] for row in rows}),
               "manifest_rows": len(rows), "audit_rows": len(audit),
               "exception_rows": len(exceptions), "strict": args.strict, "dry_run": args.dry_run,
               "published": bool(rows and not exceptions and not args.dry_run)}
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    output_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if exceptions and args.strict:
        return 2
    if not args.dry_run:
        write_tsv(args.output_manifest, MANIFEST_FIELDS, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
