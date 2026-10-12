#!/usr/bin/env python3
"""Build a Webin V2 MODIFY document adding child studies to an umbrella."""

import argparse
import csv
import os
import xml.etree.ElementTree as ET


def mk(parent, tag, text=None, **attrib):
    element = ET.SubElement(parent, tag, attrib)
    if text is not None:
        element.text = text
    return element


def read_child_accessions(path):
    accessions = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = (row.get("accession") or "").strip()
            status = (row.get("status") or "").strip().upper()
            if accession and status in {"SUCCESS", "EXISTS"} and accession.startswith("PRJ"):
                accessions.append(accession)
    return list(dict.fromkeys(accessions))


def build_xml(umbrella_accession, child_accessions, title, description):
    root = ET.Element("WEBIN")
    submission_set = mk(root, "SUBMISSION_SET")
    submission = mk(submission_set, "SUBMISSION")
    actions = mk(submission, "ACTIONS")
    action = mk(actions, "ACTION")
    mk(action, "MODIFY")

    project_set = mk(root, "PROJECT_SET")
    project = mk(project_set, "PROJECT", accession=umbrella_accession)
    mk(project, "TITLE", title)
    mk(project, "DESCRIPTION", description)
    mk(project, "UMBRELLA_PROJECT")
    related = mk(project, "RELATED_PROJECTS")
    for accession in child_accessions:
        relation = mk(related, "RELATED_PROJECT")
        mk(relation, "CHILD_PROJECT", accession=accession)
    return ET.ElementTree(root)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--umbrella-accession", required=True)
    parser.add_argument("--child-accessions", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--description", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    children = read_child_accessions(args.child_accessions)
    if not children:
        raise RuntimeError("no successful child project accessions found")
    if len(args.title) < 20 or len(args.description) < 20:
        raise RuntimeError("umbrella title and description must each contain at least 20 characters")
    os.makedirs(args.outdir, exist_ok=True)
    path = os.path.join(args.outdir, "umbrella.xml")
    build_xml(args.umbrella_accession, children, args.title, args.description).write(
        path, encoding="UTF-8", xml_declaration=True
    )
    print(path)


if __name__ == "__main__":
    main()
