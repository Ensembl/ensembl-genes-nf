#!/usr/bin/env python3
"""Prepare a safe ENA ANALYSIS MODIFY payload.

The input XML should be the current analysis XML retrieved from ENA.  The
script preserves the existing file declarations and replaces only SAMPLE_REF
elements with ENA Sample accessions resolved from the source runs.
"""

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

try:
    from .resolve_ena_metadata import (
        ENA_SAMPLE_RE,
        BIO_SAMPLE_RE,
        canonical_sample,
        fetch_run,
        unique,
    )
except ImportError:  # pragma: no cover - direct script execution
    from resolve_ena_metadata import (  # type: ignore
        ENA_SAMPLE_RE,
        BIO_SAMPLE_RE,
        canonical_sample,
        fetch_run,
        unique,
    )


def find_analysis(root):
    if root.tag == "ANALYSIS":
        return root
    analyses = root.findall(".//ANALYSIS")
    if len(analyses) != 1:
        raise RuntimeError(f"Expected exactly one ANALYSIS element, found {len(analyses)}")
    return analyses[0]


def read_current_analysis(path):
    tree = ET.parse(path)
    analysis = find_analysis(tree.getroot())
    accession = (analysis.get("accession") or "").strip()
    alias = (analysis.get("alias") or "").strip()
    if not accession and not alias:
        raise RuntimeError("Current analysis XML has neither accession nor alias")
    files = analysis.findall("FILES/FILE")
    if not files:
        raise RuntimeError("Current analysis XML has no FILE declarations")
    return tree, analysis


def source_runs(analysis):
    return unique(
        ref.get("accession", "").strip()
        for ref in analysis.findall("RUN_REF")
        if ref.get("accession")
    )


def resolve_samples(analysis, runs, explicit_samples, api_url):
    if explicit_samples and runs:
        raise RuntimeError("Use --sample-accession or --run-accession, not both")
    if explicit_samples:
        invalid = [value for value in explicit_samples if not ENA_SAMPLE_RE.match(value)]
        biosamples = [value for value in explicit_samples if BIO_SAMPLE_RE.match(value)]
        if invalid:
            raise RuntimeError(
                "--sample-accession accepts ENA Sample accessions only (ERS/SRS/DRS): "
                + ", ".join(invalid)
            )
        return unique(explicit_samples), biosamples
    if not runs:
        raise RuntimeError("No RUN_REF found; provide --run-accession or --sample-accession")

    records = [(run, fetch_run(run, api_url)) for run in runs]
    unresolved = [
        f"{run}={canonical_sample(record) or '<missing>'}"
        for run, record in records
        if not ENA_SAMPLE_RE.match(canonical_sample(record))
    ]
    if unresolved:
        raise RuntimeError(
            "Could not resolve all source runs to ENA Sample accessions: "
            + ", ".join(unresolved)
        )
    samples = unique(canonical_sample(record) for _, record in records)
    biosamples = unique(
        value
        for _, record in records
        for value in (record.get("sample_accession", ""), record.get("secondary_sample_accession", ""))
        if BIO_SAMPLE_RE.match((value or "").strip())
    )
    return samples, biosamples


def replace_sample_refs(analysis, samples, biosamples=None):
    old_refs = analysis.findall("SAMPLE_REF")
    for ref in old_refs:
        analysis.remove(ref)

    biosamples = biosamples or []
    new_refs = []
    for index, sample in enumerate(samples):
        ref = ET.Element("SAMPLE_REF", accession=sample)
        if len(samples) == len(biosamples):
            identifiers = ET.SubElement(ref, "IDENTIFIERS")
            ET.SubElement(identifiers, "PRIMARY_ID").text = sample
            ET.SubElement(identifiers, "EXTERNAL_ID", namespace="BioSample").text = biosamples[index]
        new_refs.append(ref)
    children = list(analysis)
    insertion_index = next(
        (index for index, child in enumerate(children) if child.tag in {"EXPERIMENT_REF", "RUN_REF", "ANALYSIS_TYPE"}),
        len(children),
    )
    for offset, ref in enumerate(new_refs):
        analysis.insert(insertion_index + offset, ref)


def add_attribute(analysis, tag, value):
    if not value:
        return
    attrs = analysis.find("ANALYSIS_ATTRIBUTES")
    if attrs is None:
        attrs = ET.SubElement(analysis, "ANALYSIS_ATTRIBUTES")
    for attr in attrs.findall("ANALYSIS_ATTRIBUTE"):
        if attr.findtext("TAG") == tag:
            value_element = attr.find("VALUE")
            if value_element is None:
                value_element = ET.SubElement(attr, "VALUE")
            value_element.text = value
            return
    attr = ET.SubElement(attrs, "ANALYSIS_ATTRIBUTE")
    ET.SubElement(attr, "TAG").text = tag
    ET.SubElement(attr, "VALUE").text = value


def build_submission_xml(validate=False):
    root = ET.Element("SUBMISSION")
    actions = ET.SubElement(root, "ACTIONS")
    action = ET.SubElement(actions, "ACTION")
    ET.SubElement(action, "MODIFY")
    if validate:
        action = ET.SubElement(actions, "ACTION")
        ET.SubElement(action, "VALIDATE")
    return ET.ElementTree(root)


def write_tree(tree, path):
    tree.write(path, encoding="UTF-8", xml_declaration=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-xml", required=True, help="Current ENA ANALYSIS XML")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--run-accession", action="append", default=[])
    parser.add_argument("--sample-accession", action="append", default=[])
    parser.add_argument("--api-url", default="https://www.ebi.ac.uk/ena/portal/api/filereport")
    args = parser.parse_args()

    tree, analysis = read_current_analysis(args.analysis_xml)
    runs = unique(args.run_accession or source_runs(analysis))
    samples, biosamples = resolve_samples(analysis, runs, unique(args.sample_accession), args.api_url)
    replace_sample_refs(analysis, samples, biosamples)
    add_attribute(analysis, "biosample_accessions", ",".join(biosamples))
    add_attribute(analysis, "source_metadata_origin", "ena_portal_api")

    os.makedirs(args.outdir, exist_ok=True)
    analysis_path = Path(args.outdir) / "analysis_modify.xml"
    submission_path = Path(args.outdir) / "submission_modify.xml"
    validation_path = Path(args.outdir) / "submission_validate.xml"
    write_tree(tree, analysis_path)
    write_tree(build_submission_xml(), submission_path)
    write_tree(build_submission_xml(validate=True), validation_path)

    print(analysis_path)
    print(submission_path)
    print(validation_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
