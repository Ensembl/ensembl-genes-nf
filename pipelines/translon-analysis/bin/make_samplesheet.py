#!/usr/bin/env python3
"""Create an explicit cohort samplesheet from a published riboseq tree."""
import argparse
import csv
from pathlib import Path


def discover(root):
    records = {}
    for bam in root.rglob("*.bam"):
        if bam.name.endswith("Aligned.toTranscriptome.out.bam"):
            sample = bam.name.removesuffix(".Aligned.toTranscriptome.out.bam")
            records.setdefault(sample, {})["transcriptome_bam"] = bam
        elif bam.name.endswith("Aligned.sortedByCoord.out.bam"):
            sample = bam.name.removesuffix(".Aligned.sortedByCoord.out.bam")
            records.setdefault(sample, {})["genome_bam"] = bam
    rows = []
    for sample in sorted(records):
        record = records[sample]
        # Keep the default conservative: samples are not pooled unless the
        # user edits merge_group (for example, to a cohort or study ID).
        row = {"sample_id": sample, "merge_group": sample}
        for kind in ("transcriptome_bam", "genome_bam"):
            bam = record.get(kind)
            row[kind] = str(bam.resolve()) if bam else ""
            row[kind.replace("bam", "bai")] = str(Path(f"{bam}.bai").resolve()) if bam and Path(f"{bam}.bai").exists() else ""
        fastqs = sorted(root.rglob(f"{sample}*.fastq.gz")) + sorted(root.rglob(f"{sample}*.fastq"))
        row["ribo_fastq"] = str(fastqs[0].resolve()) if fastqs else ""
        offsets = []
        for suffix in ("offsets.selected.tsv", "offsets.good.tsv", "offsets.pass.tsv"):
            offsets.extend(root.rglob(f"{sample}*.{suffix}"))
        row["offsets"] = str(sorted(offsets)[0].resolve()) if offsets else ""
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--riboseq-outdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = discover(args.riboseq_outdir)
    if not rows:
        raise SystemExit(f"No published STAR BAMs found below {args.riboseq_outdir}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample_id", "merge_group", "transcriptome_bam", "transcriptome_bai", "genome_bam", "genome_bai", "ribo_fastq", "offsets"]
    for row in rows:
        row.setdefault("ribo_fastq", "")
        row.setdefault("offsets", "")
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
