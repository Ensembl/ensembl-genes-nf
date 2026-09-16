#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("manifest")
    p.add_argument("output")
    a = p.parse_args()
    with open(a.manifest, newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        accessions = [row["run_accession"] for row in rows if row.get("run_accession")]
    Path(a.output).write_text("\n".join(accessions) + ("\n" if accessions else ""))


if __name__ == "__main__":
    main()
