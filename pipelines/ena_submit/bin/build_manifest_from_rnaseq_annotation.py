#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional


MANIFEST_COLUMNS = [
    "files_tsv",
    "study",
    "project_alias",
    "umbrella_study",
    "analysis_alias",
    "title",
    "description",
    "assembly_accession",
    "reference_fasta",
    "assembly_report",
    "reference_supplement",
    "last_geneset_update",
    "partial_release_label",
    "species",
    "taxon_id",
    "ref_seqs",
    "analysis_links",
    "analysis_attributes",
    "analysis_type",
    "omit_run_refs_in_test",
]

FILES_COLUMNS = [
    "file_path",
    "file_type",
    "remote_name",
    "run_accession",
    "sample_accession",
    "experiment_accession",
    "platform",
    "source_fastq_count",
    "source_fastq_urls",
    "source_fastq_md5s",
    "bam_sort_order",
    "bam_pg_programs",
    "alignment_software",
    "alignment_software_version",
]


def slug(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return re.sub(r"_+", "_", value).strip("_")


def partial_release_label(assembly_accession: str, last_geneset_update: str) -> str:
    return f"{assembly_accession}-Ensembl-{last_geneset_update}"


def default_species_from_rnaseq_dir(rnaseq_dir: Path) -> str:
    # Expected layout: .../<species>/<assembly_accession>/rnaseq
    try:
        return rnaseq_dir.parent.parent.name
    except Exception:
        return ""


def default_runs_csv(rnaseq_dir: Path, species: Optional[str]) -> Path:
    candidates = []
    if species:
        candidates.append(rnaseq_dir / f"{species}.csv")
    candidates.extend(sorted(p for p in rnaseq_dir.glob("*.csv") if not p.name.endswith("_gen.csv")))
    candidates.extend(sorted(rnaseq_dir.glob("*.csv")))
    seen = []
    for candidate in candidates:
        if candidate not in seen:
            seen.append(candidate)
    existing = [p for p in seen if p.exists()]
    if not existing:
        raise FileNotFoundError(f"No run metadata CSV found under {rnaseq_dir}")
    return existing[0]


def read_run_metadata(path: Path) -> OrderedDict:
    runs = OrderedDict()
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, start=1):
            if not row or all(not col.strip() for col in row):
                continue
            if len(row) < 2:
                raise ValueError(f"{path}:{line_no}: expected at least sample and run columns")
            sample = row[0].strip()
            run = row[1].strip()
            if not run:
                continue
            entry = runs.setdefault(
                run,
                {
                    "sample_accession": sample,
                    "platforms": OrderedDict(),
                    "fastq_urls": [],
                    "fastq_md5s": [],
                },
            )
            if sample and not entry["sample_accession"]:
                entry["sample_accession"] = sample
            if len(row) > 8 and row[8].strip():
                entry["platforms"][row[8].strip()] = True
            if len(row) > 10 and row[10].strip():
                entry["fastq_urls"].append(row[10].strip())
            if len(row) > 11 and row[11].strip():
                entry["fastq_md5s"].append(row[11].strip())
    return runs


def infer_file_type(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    suffix = path.suffix.lower().lstrip(".")
    if suffix in {"bam", "cram"}:
        return suffix
    raise ValueError(f"Cannot infer file type from {path}; use --file-format bam|cram")


def find_alignment_file(output_dir: Path, run: str, file_format: str, file_pattern: str) -> Optional[Path]:
    if file_format == "auto":
        for ext in ("bam", "cram"):
            candidate = output_dir / file_pattern.format(run=run, ext=ext)
            if candidate.exists():
                return candidate
        return None
    candidate = output_dir / file_pattern.format(run=run, ext=file_format)
    return candidate if candidate.exists() else None


def find_reference_files(
    assembly_dir: Path,
    assembly_accession: str,
    reference_fasta: Optional[str],
    assembly_report: Optional[str],
) -> tuple[Path, Path]:
    """Resolve the INSDC FASTA and NCBI assembly report beside the assembly."""
    if reference_fasta:
        fasta = Path(reference_fasta).expanduser().resolve()
    else:
        candidates = sorted(
            path for path in assembly_dir.glob(f"{assembly_accession}_*_genomic.fna*")
            if path.is_file() and not path.name.endswith(".fai")
        )
        if len(candidates) != 1:
            raise RuntimeError(
                f"Expected exactly one INSDC genomic FASTA under {assembly_dir}; found {candidates}. "
                "Use --reference-fasta to select one."
            )
        fasta = candidates[0].resolve()

    if assembly_report:
        report = Path(assembly_report).expanduser().resolve()
    else:
        candidates = sorted(
            path for path in assembly_dir.glob(f"{assembly_accession}_*_assembly_report.txt")
            if path.is_file()
        )
        if len(candidates) != 1:
            raise RuntimeError(
                f"Expected exactly one assembly report under {assembly_dir}; found {candidates}. "
                "Use --assembly-report to select one."
            )
        report = candidates[0].resolve()

    if not fasta.exists():
        raise FileNotFoundError(f"INSDC reference FASTA not found: {fasta}")
    if not report.exists():
        raise FileNotFoundError(f"Assembly report not found: {report}")
    return fasta, report


def assembly_report_rows(report: Path) -> list[dict[str, str]]:
    columns = None
    rows = []
    with report.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                candidate = line.lstrip("# ").split("\t")
                if "Sequence-Name" in candidate:
                    columns = candidate
                continue
            if not line.strip():
                continue
            if columns is None:
                columns = [
                    "Sequence-Name", "Sequence-Role", "Assigned-Molecule",
                    "Assigned-Molecule-Location/Type", "GenBank-Accn",
                    "Relationship", "RefSeq-Accn", "Assembly-Unit",
                    "Sequence-Length", "UCSC-style-name",
                ]
            rows.append(dict(zip(columns, line.split("\t"))) )
    return rows


def report_accessions(row: dict[str, str]) -> list[str]:
    accessions = []
    for key in ("RefSeq-Accn", "GenBank-Accn", "Sequence-Name"):
        value = (row.get(key) or "").strip()
        if value and value.lower() != "na" and value not in accessions:
            accessions.append(value)
    return accessions


def locate_supplement(assembly_dir: Path, accession: str, excluded: set[Path]) -> Optional[Path]:
    candidates = []
    for path in sorted(assembly_dir.glob(f"{accession}*")):
        if path.is_file() and path.resolve() not in excluded and not path.name.endswith(".fai"):
            candidates.append(path.resolve())
    if not candidates:
        return None
    if len(candidates) > 1:
        raise RuntimeError(f"Multiple files found for assembly supplement {accession}: {candidates}")
    return candidates[0]


def supplement_to_fasta(source: Path, output, accession: str) -> None:
    """Append a FASTA or a simple GenBank flatfile record to output."""
    opener = gzip.open if source.name.endswith(".gz") else open
    with opener(source, "rt") as handle:
        first = handle.readline()
        if first.startswith(">"):
            output.write(first)
            for line in handle:
                output.write(line)
            return
        if not first.startswith("LOCUS"):
            raise RuntimeError(f"Unsupported supplement format for {accession}: {source}")
        in_origin = False
        sequence = []
        for line in handle:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if line.startswith("//"):
                break
            if in_origin:
                sequence.append("".join(ch for ch in line if ch.isalpha()))
        sequence = "".join(sequence).upper()
        if not sequence:
            raise RuntimeError(f"No sequence found in GenBank supplement {source}")
        output.write(f">{accession}\n")
        for start in range(0, len(sequence), 80):
            output.write(sequence[start:start + 80] + "\n")


def build_reference_supplement(assembly_dir: Path, fasta: Path, report: Path, output: Path) -> list[str]:
    """Find report-listed non-nuclear records absent from the genomic FASTA."""
    genomic_names = set()
    with fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                genomic_names.add(line[1:].split()[0])

    selected = []
    excluded = {fasta.resolve(), report.resolve()}
    for row in assembly_report_rows(report):
        location_type = (row.get("Assigned-Molecule-Location/Type") or "").lower()
        role = (row.get("Sequence-Role") or "").lower()
        if "non-nuclear" not in location_type and "mitochond" not in location_type and "organel" not in role:
            continue
        accession = next((a for a in report_accessions(row) if a not in genomic_names), None)
        if not accession:
            continue
        source = locate_supplement(assembly_dir, accession, excluded)
        if source is None:
            raise RuntimeError(
                f"Assembly report lists non-nuclear sequence {accession}, but no matching record was found under {assembly_dir}"
            )
        selected.append((accession, source))
        excluded.add(source)

    with output.open("w") as handle:
        for accession, source in selected:
            supplement_to_fasta(source, handle, accession)
    return [accession for accession, _ in selected]


def parse_bam_header(path: Path) -> Dict[str, str]:
    if not shutil.which("samtools"):
        return {}
    try:
        result = subprocess.run(
            ["samtools", "view", "-H", str(path)],
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
    except Exception:
        return {}

    sort_order = ""
    programs = OrderedDict()
    star_version = ""
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if not parts:
            continue
        if parts[0] == "@HD":
            for field in parts[1:]:
                if field.startswith("SO:"):
                    sort_order = field[3:]
        elif parts[0] == "@PG":
            pg = {}
            for field in parts[1:]:
                if ":" in field:
                    key, value = field.split(":", 1)
                    pg[key] = value
            name = pg.get("PN") or pg.get("ID")
            if name:
                programs[name] = True
            if name and name.upper() == "STAR" and pg.get("VN"):
                star_version = pg["VN"]

    return {
        "bam_sort_order": sort_order,
        "bam_pg_programs": ",".join(programs.keys()),
        "alignment_software": "STAR" if star_version or "STAR" in programs else "",
        "alignment_software_version": star_version,
    }


def write_tsv(path: Path, columns: List[str], rows: List[Dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build an annotation-level ENA manifest from an Ensembl genebuild RNA-seq directory."
    )
    parser.add_argument("--rnaseq-dir", required=True, help="Path to .../<species>/<assembly>/rnaseq")
    parser.add_argument("--runs-csv", help="Headerless RNA-seq run metadata TSV; defaults to <species>.csv")
    parser.add_argument("--assembly-accession", required=True)
    parser.add_argument("--reference-fasta", help="INSDC genomic FASTA; auto-discovered beside the assembly")
    parser.add_argument("--assembly-report", help="NCBI assembly report; auto-discovered beside the assembly")
    parser.add_argument("--last-geneset-update", required=True, help="YYYY-MM value from genome metadata")
    parser.add_argument("--species", help="Production species name; defaults from directory layout")
    parser.add_argument("--taxon-id", default="")
    parser.add_argument("--study", default="", help="Existing child study accession/alias. Usually blank when auto-creating.")
    parser.add_argument("--project-alias", help="Child project alias. Defaults to prj_<partial release label>.")
    parser.add_argument("--umbrella-study", default="", help="Umbrella study accession/alias, recorded as metadata.")
    parser.add_argument("--file-format", choices=["auto", "bam", "cram"], default="auto")
    parser.add_argument(
        "--file-pattern",
        default="{run}_Aligned.sortedByCoord.out.{ext}",
        help="Pattern under output/ with {run} and {ext} placeholders.",
    )
    parser.add_argument("--output-subdir", default="output")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--analysis-alias", help="Stable ENA analysis alias.")
    parser.add_argument("--title", help="Analysis title.")
    parser.add_argument("--description", help="Analysis description.")
    parser.add_argument("--analysis-links", default="")
    parser.add_argument("--extra-attribute", action="append", default=[], help="Extra attr key=value; repeatable.")
    parser.add_argument("--allow-missing-files", action="store_true")
    parser.add_argument("--no-bam-header", action="store_true", help="Do not inspect BAM/CRAM headers with samtools.")
    args = parser.parse_args()

    rnaseq_dir = Path(args.rnaseq_dir).resolve()
    assembly_dir = rnaseq_dir.parent
    reference_fasta, assembly_report = find_reference_files(
        assembly_dir, args.assembly_accession, args.reference_fasta, args.assembly_report
    )
    output_dir = rnaseq_dir / args.output_subdir
    if not output_dir.exists():
        raise FileNotFoundError(f"Alignment output directory not found: {output_dir}")

    species = args.species or default_species_from_rnaseq_dir(rnaseq_dir)
    runs_csv = Path(args.runs_csv).resolve() if args.runs_csv else default_runs_csv(rnaseq_dir, species)
    runs = read_run_metadata(runs_csv)
    if not runs:
        raise RuntimeError(f"No runs found in {runs_csv}")

    release_label = partial_release_label(args.assembly_accession, args.last_geneset_update)
    analysis_alias = args.analysis_alias or slug(
        f"rnaseq_alignment_evidence_{args.assembly_accession}_Ensembl_{args.last_geneset_update}"
    )
    project_alias = args.project_alias or slug(f"prj_{release_label}")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    reference_supplement = outdir / "reference_supplement.fasta"
    supplement_accessions = build_reference_supplement(
        assembly_dir, reference_fasta, assembly_report, reference_supplement
    )
    files_tsv = outdir / "files.tsv"
    manifest_tsv = outdir / "manifest.tsv"
    missing_tsv = outdir / "missing_files.tsv"
    summary_tsv = outdir / "summary.tsv"

    file_rows = []
    missing_rows = []
    for run, info in runs.items():
        alignment_file = find_alignment_file(output_dir, run, args.file_format, args.file_pattern)
        if alignment_file is None:
            missing_rows.append({"run_accession": run, "sample_accession": info["sample_accession"]})
            continue
        file_type = infer_file_type(alignment_file, args.file_format)
        remote_name = f"{run}_{release_label}.{file_type}"
        header_meta = {} if args.no_bam_header else parse_bam_header(alignment_file)
        file_rows.append(
            {
                "file_path": str(alignment_file),
                "file_type": file_type,
                "remote_name": remote_name,
                "run_accession": run,
                "sample_accession": info["sample_accession"],
                "experiment_accession": "",
                "platform": ",".join(info["platforms"].keys()),
                "source_fastq_count": str(len(info["fastq_urls"])),
                "source_fastq_urls": ",".join(info["fastq_urls"]),
                "source_fastq_md5s": ",".join(info["fastq_md5s"]),
                **header_meta,
            }
        )

    if missing_rows and not args.allow_missing_files:
        write_tsv(missing_tsv, ["run_accession", "sample_accession"], missing_rows)
        raise RuntimeError(f"{len(missing_rows)} runs had no alignment file; see {missing_tsv}")
    if not file_rows:
        raise RuntimeError("No alignment files found for submission")

    software = sorted({row["alignment_software"] for row in file_rows if row.get("alignment_software")})
    software_versions = sorted(
        {row["alignment_software_version"] for row in file_rows if row.get("alignment_software_version")}
    )
    platforms = sorted({p for row in file_rows for p in row.get("platform", "").split(",") if p})

    attributes = OrderedDict(
        [
            ("pipeline", "ensembl-genebuild-rnaseq"),
            ("evidence_type", "RNA-seq alignments"),
            ("processing_level", "processed alignments assessed for annotation evidence, no annotation filtering decisions"),
            ("assembly_accession", args.assembly_accession),
            ("ensembl_partial_release", release_label),
            ("last_geneset_update", args.last_geneset_update),
            ("species", species),
            ("source_runs_count", str(len(file_rows))),
            ("alignment_files_count", str(len(file_rows))),
            ("source_metadata_origin", "rnaseq_csv"),
        ]
    )
    if args.taxon_id:
        attributes["taxon_id"] = args.taxon_id
    if args.umbrella_study:
        attributes["umbrella_study"] = args.umbrella_study
    if platforms:
        attributes["platform"] = ",".join(platforms)
    if software:
        attributes["alignment_software"] = ",".join(software)
    if software_versions:
        attributes["alignment_software_version"] = ",".join(software_versions)
    for item in args.extra_attribute:
        if "=" not in item:
            raise ValueError(f"--extra-attribute must be key=value: {item}")
        key, value = item.split("=", 1)
        attributes[key.strip()] = value.strip()

    title = args.title or f"RNA-seq alignment evidence for {release_label}"
    description = args.description or (
        f"Processed RNA-seq alignments assessed as evidence for Ensembl annotation of "
        f"{args.assembly_accession}; last geneset update {args.last_geneset_update}."
    )
    manifest_row = {
        "files_tsv": str(files_tsv),
        "study": args.study,
        "project_alias": project_alias,
        "umbrella_study": args.umbrella_study,
        "analysis_alias": analysis_alias,
        "title": title,
        "description": description,
        "assembly_accession": args.assembly_accession,
        "reference_fasta": str(reference_fasta),
        "assembly_report": str(assembly_report),
        "reference_supplement": str(reference_supplement),
        "last_geneset_update": args.last_geneset_update,
        "partial_release_label": release_label,
        "species": species,
        "taxon_id": args.taxon_id,
        "ref_seqs": ",".join(supplement_accessions),
        "analysis_links": args.analysis_links,
        "analysis_attributes": "; ".join(f"attr_{key}={value}" for key, value in attributes.items() if value),
        "analysis_type": "REFERENCE_ALIGNMENT",
        "omit_run_refs_in_test": "true",
    }

    write_tsv(files_tsv, FILES_COLUMNS, file_rows)
    write_tsv(manifest_tsv, MANIFEST_COLUMNS, [manifest_row])
    write_tsv(missing_tsv, ["run_accession", "sample_accession"], missing_rows)
    write_tsv(
        summary_tsv,
        ["metric", "value"],
        [
            {"metric": "runs_in_csv", "value": str(len(runs))},
            {"metric": "alignment_files", "value": str(len(file_rows))},
            {"metric": "missing_alignment_files", "value": str(len(missing_rows))},
            {"metric": "partial_release_label", "value": release_label},
            {"metric": "analysis_alias", "value": analysis_alias},
            {"metric": "project_alias", "value": project_alias},
        ],
    )
    print(manifest_tsv)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
