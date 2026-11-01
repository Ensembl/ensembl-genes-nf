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
EXCEPTION_FIELDS = ["source_id", "source_type", "entity_id", "gca_accession", "exception_class", "path", "detail"]
GCA_RE = re.compile(r"GCA_[0-9]+\.[0-9]+")
RUN_RE = re.compile(r"(?:ERR|SRR|DRR)[0-9]+")
BAM_SUFFIX = "_Aligned.sortedByCoord.out.bam"
GTF_SUFFIXES = (
    "_Aligned.sortedByCoord.out.gtf",
    ".scallop.gtf",
    ".stringtie.gtf",
    ".gtf",
)
GDM_IGNORED_ROOTS = {"cache", "ncbi_dataset", "nfs", "pipeline_info"}


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=("anno", "gdm", "main"))
    parser.add_argument("--base-dir", required=True, type=Path)
    parser.add_argument("--entity", help="Optional run ID or directory below base-dir")
    parser.add_argument("--annotation-dir", type=Path, help="Override the mode's annotation directory")
    parser.add_argument("--gca", help="GCA accession; inferred from the base path when omitted")
    parser.add_argument("--assembly-release", help="Assembly release; defaults to --gca")
    parser.add_argument("--gdm-duplicate-policy", choices=("error", "newest"), default="newest",
                        help="How to resolve repeated GDM run BAMs (default: newest unique mtime within assembly)")
    parser.add_argument("--gene-model-source", choices=("manifest", "scallop"), default="manifest",
                        help="Use manifest GTFs or create run-level GTFs with Scallop")
    parser.add_argument(
        "--track-types",
        help="Comma-delimited tracks; defaults to all for main/anno, or coverage+splice_junction for gdm (all three with --gene-model-source scallop)",
    )
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--output-manifest", type=Path)
    output.add_argument("--output-dir", type=Path,
                        help="Write one manifest and diagnostics per GCA under this directory")
    parser.add_argument("--output-audit", type=Path)
    parser.add_argument("--output-exceptions", type=Path)
    parser.add_argument("--output-summary", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-strict", dest="strict", action="store_false")
    parser.set_defaults(strict=True)
    return parser.parse_args(argv)


def normalise_tracks(raw: Optional[str], mode: str, gene_model_source: str = "manifest") -> List[str]:
    if raw is None or not raw.strip():
        raw = ("coverage,gene_model,splice_junction"
               if mode == "gdm" and gene_model_source == "scallop"
               else ("coverage,splice_junction" if mode == "gdm" else "all"))
    selected = [part.strip() for part in raw.split(",") if part.strip()]
    if not selected:
        selected = (
            (["coverage", "gene_model", "splice_junction"]
             if mode == "gdm" and gene_model_source == "scallop"
             else ["coverage", "splice_junction"])
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
    tracks = normalise_tracks(args.track_types, args.mode, args.gene_model_source)
    required = {"coverage": bam,
                "gene_model": True if args.gene_model_source == "scallop" else gtf,
                "splice_junction": sj}
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


def gdm_inventory_files(base: Path) -> Tuple[List[Path], List[Dict[str, str]]]:
    """Find GDM report/readme inventories while retaining permission failures."""
    files: List[Path] = []
    exceptions: List[Dict[str, str]] = []

    def onerror(error: OSError) -> None:
        exceptions.append({
            "source_id": base.name,
            "source_type": "gdm_inventory",
            "entity_id": "",
            "gca_accession": "",
            "exception_class": "unreadable_directory",
            "path": str(error.filename or base),
            "detail": str(error),
        })

    for root, directories, names in os.walk(base, onerror=onerror):
        directories[:] = sorted(name for name in directories if name not in GDM_IGNORED_ROOTS)
        for name in sorted(names):
            if name in {"readme.csv", "report.csv"}:
                files.append(Path(root) / name)
    return sorted(files), exceptions


def gdm_species_root(inventory: Path, base: Path) -> Path:
    for parent in (inventory.parent, *inventory.parents):
        if re.search(r"(?i)_gca[0-9]+v[0-9]+$", parent.name):
            return parent
        if parent == base:
            break
    return inventory.parent


def resolve_gdm_inventory_path(raw: str, inventory: Path, base: Path) -> Optional[Path]:
    if not raw or not raw.strip():
        return None
    candidate = Path(raw)
    if candidate.is_absolute():
        return candidate if candidate.is_file() else None
    species_root = gdm_species_root(inventory, base)
    candidates = [species_root / candidate, inventory.parent / candidate, base / candidate]
    for path in candidates:
        if path.is_file():
            return path
    return None


def gdm_exception(args: argparse.Namespace, inventory: Path, entity_id: str,
                  gca: str, exception_class: str, detail: str, path: Optional[Path] = None) -> Dict[str, str]:
    return {
        "source_id": inventory.parent.name,
        "source_type": "gdm_inventory",
        "entity_id": entity_id,
        "gca_accession": gca,
        "exception_class": exception_class,
        "path": str(path or inventory),
        "detail": detail,
    }


def discover_gdm_inventory(args: argparse.Namespace, base: Path, inventory_files: Sequence[Path],
                           scan_exceptions: Sequence[Dict[str, str]]) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], List[Dict[str, str]]]:
    rows: List[Dict[str, str]] = []
    audit: List[Dict[str, str]] = []
    exceptions: List[Dict[str, str]] = list(scan_exceptions)
    candidates: List[Dict[str, object]] = []

    for inventory in inventory_files:
        try:
            with inventory.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                required = {"assembly_accession", "run_accession", "bam_file", "splice_junction_file"}
                missing = sorted(required - set(reader.fieldnames or []))
                if missing:
                    exceptions.append(gdm_exception(
                        args, inventory, "", "", "invalid_inventory",
                        "missing columns: " + ", ".join(missing)))
                    continue
                for line_number, record in enumerate(reader, start=2):
                    run_id = (record.get("run_accession") or "").strip()
                    gca = (record.get("assembly_accession") or "").strip()
                    if not RUN_RE.fullmatch(run_id):
                        exceptions.append(gdm_exception(
                            args, inventory, run_id, gca, "invalid_run_accession",
                            f"line {line_number} has an invalid run_accession", inventory))
                        continue
                    if not GCA_RE.fullmatch(gca):
                        exceptions.append(gdm_exception(
                            args, inventory, run_id, gca, "assembly_mismatch",
                            f"line {line_number} has an invalid assembly_accession", inventory))
                        continue
                    if args.gca and args.gca != gca:
                        continue
                    if args.entity and not (absolute(args.base_dir) / args.entity).is_dir() and args.entity != run_id:
                        continue
                    bam = resolve_gdm_inventory_path(record.get("bam_file", ""), inventory, base)
                    sj = resolve_gdm_inventory_path(record.get("splice_junction_file", ""), inventory, base)
                    if bam is None:
                        exceptions.append(gdm_exception(args, inventory, run_id, gca, "missing_bam",
                                                         f"line {line_number} bam_file was not found"))
                        continue
                    if sj is None:
                        exceptions.append(gdm_exception(args, inventory, run_id, gca, "missing_splice_junction",
                                                         f"line {line_number} splice_junction_file was not found", bam))
                        continue
                    bai = Path(str(bam) + ".bai") if Path(str(bam) + ".bai").is_file() else None
                    if bai is None and Path(str(bam) + ".csi").is_file():
                        bai = Path(str(bam) + ".csi")
                    candidates.append({
                        "inventory": inventory, "run_id": run_id, "gca": gca,
                        "bam": bam, "sj": sj, "bai": bai, "line_number": line_number,
                    })
        except OSError as exc:
            exceptions.append(gdm_exception(args, inventory, "", "", "unreadable_inventory", str(exc)))

    grouped: Dict[Tuple[str, str], List[Dict[str, object]]] = {}
    for candidate in candidates:
        grouped.setdefault((str(candidate["gca"]), str(candidate["run_id"])), []).append(candidate)

    for (gca, run_id), group in sorted(grouped.items()):
        unique: Dict[str, Dict[str, object]] = {str(item["bam"]): item for item in group}
        warning = ""
        if len(unique) > 1:
            if args.gdm_duplicate_policy == "error":
                exceptions.append(gdm_exception(
                    args, group[0]["inventory"], run_id, gca, "duplicate_entity",
                    f"{len(unique)} distinct BAM paths resolve to this run", group[0]["bam"]))
                continue
            ordered = sorted(unique.values(), key=lambda item: (Path(str(item["bam"])).stat().st_mtime_ns, str(item["bam"])))
            newest = ordered[-1]
            if len(ordered) > 1 and Path(str(ordered[-2]["bam"])).stat().st_mtime_ns == Path(str(newest["bam"])).stat().st_mtime_ns:
                exceptions.append(gdm_exception(
                    args, newest["inventory"], run_id, gca, "duplicate_entity",
                    "duplicate BAM candidates have tied modification times", newest["bam"]))
                continue
            selected = newest
            warning = f"selected newest of {len(unique)} candidate BAM paths by mtime"
        else:
            selected = next(iter(unique.values()))
            if len(group) > 1:
                warning = f"collapsed {len(group)} duplicate inventory rows pointing to the same BAM"

        inventory = selected["inventory"]
        try:
            row = make_row(args, selected["bam"], None, selected["sj"], selected["bai"], run_id, gca)
        except ValueError as exc:
            detail = str(exc)
            exception_class = "missing_" + detail[len("missing required "):].split(",", 1)[0] \
                if detail.startswith("missing required ") else "invalid_inventory"
            exceptions.append(gdm_exception(args, inventory, run_id, gca, exception_class, detail, selected["bam"]))
            continue
        rows.append(row)
        audit.append({
            "source_id": inventory.parent.name,
            "source_type": "gdm_inventory",
            "entity_id": run_id,
            "entity_type": "run",
            "gca_accession": gca,
            "source_path": str(absolute(inventory)),
            "bam": row["bam"], "bai": row["bai"], "gtf": row["gtf"], "sj_out_tab": row["sj_out_tab"],
            "matching_rule": inventory.name,
            "status": "resolved", "warning": warning,
        })

    rows.sort(key=lambda row: (row["gca_accession"], row["entity_id"]))
    audit.sort(key=lambda row: (row["gca_accession"], row["entity_id"]))
    if not rows and not exceptions:
        exceptions.append(gdm_exception(args, base / "report.csv", "", "", "unsupported_layout",
                                        "GDM inventories contained no resolvable run rows", base))
    return rows, audit, exceptions


def discover(args: argparse.Namespace) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], List[Dict[str, str]]]:
    base, annotation = mode_roots(args)
    if args.mode == "gdm":
        inventory_files, scan_exceptions = gdm_inventory_files(base)
        if inventory_files:
            return discover_gdm_inventory(args, base, inventory_files, scan_exceptions)
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


def write_single_output(args: argparse.Namespace, rows: Sequence[Dict[str, str]],
                        audit: Sequence[Dict[str, str]], exceptions: Sequence[Dict[str, str]]) -> None:
    assert args.output_manifest is not None
    output_manifest = args.output_manifest
    output_audit = args.output_audit or output_manifest.with_name("manifest_audit.tsv")
    output_exceptions = args.output_exceptions or output_manifest.with_name("manifest_exceptions.tsv")
    output_summary = args.output_summary or output_manifest.with_name("manifest_summary.json")
    write_tsv(output_audit, AUDIT_FIELDS, audit)
    write_tsv(output_exceptions, EXCEPTION_FIELDS, exceptions)
    summary = {"mode": args.mode, "base_dir": str(absolute(args.base_dir)),
               "gca_accessions": sorted({row["gca_accession"] for row in rows}),
               "manifest_rows": len(rows), "audit_rows": len(audit),
               "exception_rows": len(exceptions), "strict": args.strict, "dry_run": args.dry_run,
               "gene_model_source": args.gene_model_source,
               "published": bool(rows and not exceptions and not args.dry_run)}
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    output_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if not args.dry_run and not (args.strict and exceptions):
        write_tsv(output_manifest, MANIFEST_FIELDS, rows)


def write_split_outputs(args: argparse.Namespace, rows: Sequence[Dict[str, str]],
                        audit: Sequence[Dict[str, str]], exceptions: Sequence[Dict[str, str]]) -> None:
    assert args.output_dir is not None
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    global_exceptions = output_dir / "manifest_exceptions.tsv"
    write_tsv(global_exceptions, EXCEPTION_FIELDS, exceptions)
    can_publish = not args.dry_run and not (args.strict and exceptions)
    index_rows = []
    gca_values = sorted({row["gca_accession"] for row in rows})
    for gca in gca_values:
        gca_dir = output_dir / gca
        gca_rows = [row for row in rows if row["gca_accession"] == gca]
        gca_audit = [item for item in audit if item["gca_accession"] == gca]
        gca_exceptions = [item for item in exceptions if item.get("gca_accession") == gca]
        manifest = gca_dir / "run_manifest.tsv"
        audit_path = gca_dir / "manifest_audit.tsv"
        exceptions_path = gca_dir / "manifest_exceptions.tsv"
        summary_path = gca_dir / "manifest_summary.json"
        write_tsv(audit_path, AUDIT_FIELDS, gca_audit)
        write_tsv(exceptions_path, EXCEPTION_FIELDS, gca_exceptions)
        summary = {"mode": args.mode, "base_dir": str(absolute(args.base_dir)),
                   "gca_accessions": [gca], "manifest_rows": len(gca_rows),
                   "audit_rows": len(gca_audit), "exception_rows": len(gca_exceptions),
                   "strict": args.strict, "dry_run": args.dry_run,
                   "gene_model_source": args.gene_model_source,
                   "published": bool(gca_rows and can_publish and not gca_exceptions)}
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        if can_publish:
            write_tsv(manifest, MANIFEST_FIELDS, gca_rows)
        index_rows.append({
            "gca_accession": gca,
            "manifest": str(manifest),
            "audit": str(audit_path),
            "exceptions": str(exceptions_path),
            "manifest_rows": str(len(gca_rows)),
            "exception_rows": str(len(gca_exceptions)),
        })
    write_tsv(output_dir / "gdm_manifest_index.tsv",
              ["gca_accession", "manifest", "audit", "exceptions", "manifest_rows", "exception_rows"],
              index_rows)
    summary = {"mode": args.mode, "base_dir": str(absolute(args.base_dir)),
               "gca_accessions": gca_values, "manifest_rows": len(rows),
               "audit_rows": len(audit), "exception_rows": len(exceptions),
               "strict": args.strict, "dry_run": args.dry_run,
               "gene_model_source": args.gene_model_source,
               "published": bool(rows and can_publish and not exceptions),
               "output_dir": str(absolute(output_dir))}
    (output_dir / "manifest_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        normalise_tracks(args.track_types, args.mode, args.gene_model_source)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    rows, audit, exceptions = discover(args)
    if args.output_dir:
        write_split_outputs(args, rows, audit, exceptions)
    else:
        write_single_output(args, rows, audit, exceptions)
    if exceptions and args.strict:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
