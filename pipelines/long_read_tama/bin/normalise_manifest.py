#!/usr/bin/env python3
"""Validate the whitespace-delimited long-read manifest.

The source has free-text fields on both sides of the fixed run metadata, so
the parser anchors on the source/platform fields and parses URL and MD5 from
the right-hand edge.
"""
import re
import sys
import csv
from pathlib import Path


SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*\.(?:fastq|fq|bam|bai|pbi)(?:\.gz)?$", re.I)
MD5 = re.compile(r"^[0-9a-fA-F]{32}$")
URL = re.compile(r"^(?:https?|ftp)://")
INTEGER = re.compile(r"^-?\d+$")


def classify_read_type(filename, platform):
    lower_name = filename.lower()
    if "_subreads" in lower_name:
        return "pacbio_subread_hint", "weak filename hint only; requires artifact and header evidence"
    if "_ccs" in lower_name or "hifi" in lower_name:
        return "pacbio_ccs_hint", "weak filename hint only; requires artifact and header evidence"
    return "long_read_unspecified", f"platform is {platform}; filename is not used as a final read classification"


def fail(message):
    raise ValueError(message)


def parse(line, number):
    tokens = line.split()
    if len(tokens) < 12:
        fail(f"line {number}: expected at least 12 whitespace-delimited fields")
    checksum, url = tokens[-1], tokens[-2]
    if not MD5.fullmatch(checksum):
        fail(f"line {number}: invalid MD5 checksum: {checksum}")
    if not URL.match(url):
        if url.startswith("ftp."):
            url = "ftp://" + url
        else:
            fail(f"line {number}: invalid download URL: {url}")
    try:
        source_index = tokens.index("ENA")
        # The historical platform is provenance only.  Read representation is
        # resolved by the inventory/classification phase, never here.
        declared_platform = tokens[source_index + 1]
        run = tokens[source_index - 6]
        filename = tokens[source_index - 4]
        numeric_positions = [source_index - 5, source_index - 3, source_index - 2, source_index - 1]
        if not all(INTEGER.fullmatch(tokens[position]) for position in numeric_positions) or not re.fullmatch(r"(?:SRR|ERR|DRR)[0-9]+", run):
            raise ValueError
    except (ValueError, IndexError):
        fail(f"line {number}: could not locate fixed ENA/platform metadata")
    if not SAFE_NAME.fullmatch(filename):
        fail(f"line {number}: unsafe or unexpected FASTQ filename: {filename}")
    tissue = " ".join(tokens[: source_index - 6]).strip()
    description = " ".join(tokens[source_index + 2 : -2]).strip()
    if not tissue or not description:
        fail(f"line {number}: tissue and description must be non-empty")
    read_type, evidence = classify_read_type(filename, declared_platform)
    return [run, tissue, filename, url, checksum.lower(), "ENA", declared_platform, description, read_type, evidence]


def main(source, destination, report):
    rows = []
    seen = set()
    errors = []
    source_lines = Path(source).read_text().splitlines()
    first_data = next((line for line in source_lines if line.strip() and not line.lstrip().startswith("#")), "")
    # Accept the versioned tabular contract in either TSV or conventional CSV
    # form.  The workflow continues to emit TSV internally, so this is an
    # input convenience rather than a second downstream channel contract.
    if "run_accession" in first_data and ("\t" in first_data or "," in first_data):
        delimiter = "\t" if "\t" in first_data else ","
        with Path(source).open(newline="") as handle:
            tab_rows = csv.DictReader(handle, delimiter=delimiter)
            required = {"run_accession", "tissue", "description", "url", "md5", "platform"}
            if not required.issubset(tab_rows.fieldnames or set()):
                errors.append("versioned TSV manifest is missing required fields")
            else:
                for number, item in enumerate(tab_rows, 2):
                    try:
                        run, url, checksum = item["run_accession"], item["url"], item["md5"]
                        if url.startswith("ftp."): url = "ftp://" + url
                        if not re.fullmatch(r"(?:SRR|ERR|DRR)[0-9]+", run) or not MD5.fullmatch(checksum) or not URL.match(url):
                            fail(f"line {number}: invalid accession, URL, or MD5")
                        filename = item.get("filename") or Path(url.split("?", 1)[0]).name
                        if not SAFE_NAME.fullmatch(filename): fail(f"line {number}: unsafe FASTQ filename: {filename}")
                        if not item["tissue"].strip() or not item["description"].strip(): fail(f"line {number}: tissue and description must be non-empty")
                        if run in seen: fail(f"line {number}: duplicate run accession: {run}")
                        seen.add(run)
                        read_type, evidence = classify_read_type(filename, item["platform"])
                        rows.append([run, item["tissue"].strip(), filename, url, checksum.lower(), item.get("source", "candidate"), item["platform"], item["description"].strip(), read_type, evidence])
                    except (KeyError, ValueError) as exc:
                        errors.append(str(exc))
        source_lines = []
    for number, raw in enumerate(source_lines, 1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        try:
            row = parse(raw, number)
            if row[0] in seen:
                fail(f"line {number}: duplicate run accession: {row[0]}")
            seen.add(row[0])
            rows.append(row)
        except ValueError as exc:
            errors.append(str(exc))
    Path(report).write_text("status\tdetail\n" + ("ok\tvalidated\n" if not errors else "error\t" + " | ".join(errors) + "\n"))
    if errors:
        raise SystemExit("\n".join(errors))
    with Path(destination).open("w") as handle:
        handle.write("run_accession\ttissue\tfilename\turl\tmd5\tsource\tplatform\tdescription\tread_type\tclassification_evidence\n")
        handle.writelines("\t".join(row) + "\n" for row in rows)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("usage: normalise_manifest.py INPUT OUTPUT REPORT")
    main(*sys.argv[1:])
