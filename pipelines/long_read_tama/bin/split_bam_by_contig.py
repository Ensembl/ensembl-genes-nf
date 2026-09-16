#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("bam")
    p.add_argument("workload")
    p.add_argument("accession")
    p.add_argument("shard_dir")
    p.add_argument("manifest")
    a = p.parse_args()
    out_dir = Path(a.shard_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(a.workload, newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        records = list(rows)
    if not records:
        raise SystemExit("No contig shards were materialised")
    manifest_rows = []
    for row in records:
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", row["contig"])
        bam = out_dir / f"{a.accession}.{safe}.bam"
        subprocess.run(["samtools", "view", "-b", "-o", str(bam), a.bam, row["contig"]], check=True)
        subprocess.run(["samtools", "index", str(bam)], check=True)
        manifest_rows.append([row["contig"], row["reference_bases"], row["mapped_reads"],
                              row["resource_class"], bam.name])
    with open(a.manifest, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "reference_bases", "mapped_reads", "resource_class", "path"])
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
