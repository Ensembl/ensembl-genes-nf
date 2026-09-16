#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("probe")
    p.add_argument("validation")
    p.add_argument("run_accession")
    p.add_argument("classification")
    p.add_argument("output")
    a = p.parse_args()
    probe = json.loads(Path(a.probe).read_text())
    validation = json.loads(Path(a.validation).read_text())
    fields = ["run_accession", "classification", "header_representation", "records_sampled",
              "distinct_ids", "distinct_molecules", "malformed", "status"]
    values = [a.run_accession, a.classification or "UNCLASSIFIED",
              probe.get("header_representation", "UNKNOWN"), probe.get("records_sampled", 0),
              validation.get("distinct_ids", 0), validation.get("distinct_molecules", 0),
              probe.get("malformed_count", 0), validation.get("status", "VALIDATED")]
    Path(a.output).write_text("\t".join(fields) + "\n" + "\t".join(map(str, values)) + "\n")


if __name__ == "__main__":
    main()
