#!/usr/bin/env python3
"""Build a merge-only manifest from completed wave output directories."""
import argparse
import csv
from pathlib import Path


FIELDS = [
    "species_id",
    "assembly_id",
    "assembly_accession",
    "run_id",
    "sample_id",
    "bam_path",
    "gtf_path",
    "assembly_status",
]


def has_transcript_records(gtf: Path) -> bool:
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] == "transcript":
                return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        action="append",
        required=True,
        help="Completed wave output directory; repeat once per wave",
    )
    parser.add_argument("--output", required=True, help="Path for the combined merge manifest")
    args = parser.parse_args()

    records = []
    seen = set()
    for outdir_arg in args.outdir:
        outdir = Path(outdir_arg).expanduser().resolve(strict=True)
        per_run = outdir / "per_run"
        if not per_run.is_dir():
            raise SystemExit(f"No per_run directory in wave output: {outdir}")
        for bam in sorted(per_run.glob("*/*/*/alignment/aligned.sorted.bam")):
            species_id, assembly_id, run_id = bam.relative_to(per_run).parts[:3]
            key = (species_id, assembly_id, run_id)
            if key in seen:
                raise SystemExit(f"Duplicate run output across wave directories: {key}")
            seen.add(key)
            assembly_dir = per_run / species_id / assembly_id / run_id / "stringtie3" / "assembly"
            gtf = assembly_dir / f"{run_id}.assembled.gtf"
            status_file = assembly_dir / "assembly_status.txt"
            if not gtf.is_file():
                status = "ASSEMBLY_MISSING"
                gtf_path = ""
            else:
                status = status_file.read_text().strip() if status_file.is_file() else ""
                if status not in {"TRANSCRIPTS_FOUND", "NO_TRANSCRIPTS"}:
                    status = "TRANSCRIPTS_FOUND" if has_transcript_records(gtf) else "NO_TRANSCRIPTS"
                gtf_path = str(gtf.resolve())
            records.append(
                {
                    "species_id": species_id,
                    "assembly_id": assembly_id,
                    "assembly_accession": assembly_id,
                    "run_id": run_id,
                    "sample_id": run_id,
                    "bam_path": str(bam.resolve()),
                    "gtf_path": gtf_path,
                    "assembly_status": status,
                }
            )

    output = Path(args.output).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    counts = {status: sum(row["assembly_status"] == status for row in records) for status in ("TRANSCRIPTS_FOUND", "NO_TRANSCRIPTS", "ASSEMBLY_MISSING")}
    print(f"Wrote {len(records)} runs to {output.resolve()}")
    print("Assembly status counts: " + ", ".join(f"{status}={count}" for status, count in counts.items()))


if __name__ == "__main__":
    main()
