#!/usr/bin/env python3
"""Left join Diamond best hits to the complete combined-model manifest."""

import csv
import argparse
import sys


CLASSIFICATIONS = {
    "long_read": [("1", 95, 90), ("2", 90, 80), ("3", 90, 60),
                   ("4", 90, 40), ("5", 80, 20), ("6", 60, 20), ("7", 0, 0)],
    "standard_old": [("1", 95, 95), ("2", 95, 80), ("3", 90, 80),
                      ("4", 80, 60), ("5", 70, 60), ("6", 50, 50),
                      ("7", 50, 25), ("8", 0, 0)],
    "standard": [("1", 95, 90), ("2", 90, 80), ("3", 90, 60),
                  ("4", 90, 40), ("5", 80, 20), ("6", 60, 20), ("7", 0, 0)],
    "gifts": [("1", 100, 100), ("2", 97, 97), ("3", 95, 95),
              ("4", 90, 90), ("5", 80, 90), ("6", 80, 80),
              ("7", 70, 70), ("8", 50, 50), ("9", 0, 0)],
}


def classify(coverage, identity, classification_type):
    for label, minimum_coverage, minimum_identity in CLASSIFICATIONS[classification_type]:
        if coverage >= minimum_coverage and identity >= minimum_identity:
            return label
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest")
    parser.add_argument("hits")
    parser.add_argument("report")
    parser.add_argument("summary")
    parser.add_argument("--classification-type", choices=CLASSIFICATIONS, default="long_read")
    args = parser.parse_args()

    manifest_path = args.manifest
    hits_path = args.hits
    report_path = args.report
    summary_path = args.summary
    hits = {}
    with open(hits_path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 13 and row[0] not in hits:
                hits[row[0]] = row

    fields = ["model_id", "query_peptide_id", "best_subject_id", "pident", "aligned_length",
              "query_coverage", "evalue", "bitscore", "legacy_classification",
              "legacy_biotype_suffix", "status"]
    statuses = {}
    with open(manifest_path) as manifest, open(report_path, "w", newline="") as report:
        reader = csv.DictReader(manifest, delimiter="\t")
        writer = csv.DictWriter(report, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            model_id = row["model_id"]
            peptide_id = row["peptide_id"]
            hit = hits.get(peptide_id) if peptide_id else None
            if hit:
                identity = float(hit[2])
                coverage = float(hit[7])
                classification = classify(coverage, identity, args.classification_type)
                status = "CLASSIFIED"
                output = [model_id, hit[0], hit[1], hit[2], hit[3], hit[7], hit[11], hit[12],
                          classification, f"_{classification}", status]
            elif row["orf_status"] == "NO_ATG":
                status = "NO_PEPTIDE"
                output = [model_id, "", "", "", "", "", "", "", "", "", status]
            else:
                status = "NO_PROTEIN_HIT"
                output = [model_id, peptide_id, "", "", "", "", "", "", "", "", status]
            writer.writerow(dict(zip(fields, output)))
            statuses[status] = statuses.get(status, 0) + 1

    with open(summary_path, "w") as summary:
        summary.write("status\tmodels\n")
        for status in sorted(statuses):
            summary.write(f"{status}\t{statuses[status]}\n")


if __name__ == "__main__":
    main()
