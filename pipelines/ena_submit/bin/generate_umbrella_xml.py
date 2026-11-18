#!/usr/bin/env python3
"""Generate ENA umbrella-project ADD or MODIFY Webin XML."""

import argparse
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET


def mk(parent, tag, text=None, **attributes):
    element = ET.SubElement(parent, tag, attributes)
    if text is not None:
        element.text = text
    return element


def child_accessions(paths):
    found = []
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                accession = (row.get("accession") or "").strip()
                status = (row.get("status") or "").strip().upper()
                if status in {"SUCCESS", "EXISTS"} and re.fullmatch(r"PRJ[A-Z]+\d+", accession):
                    if accession not in found:
                        found.append(accession)
    if not found:
        raise RuntimeError("No successful child PRJ accession found in project poll output")
    return found


def validate_analysis(paths):
    seen = False
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                status = (row.get("status") or "").strip().upper()
                if status in {"SUCCESS", "EXISTS"}:
                    seen = True
                elif status:
                    raise RuntimeError(f"Analysis submission did not succeed: {status}")
    if not seen:
        raise RuntimeError("No successful analysis accession found in analysis poll output")


def build_project(alias, title, description, children, accession=None):
    attributes = {"alias": alias}
    if accession:
        attributes["accession"] = accession
    project_set = ET.Element("PROJECT_SET")
    project = mk(project_set, "PROJECT", **attributes)
    mk(project, "TITLE", title)
    mk(project, "DESCRIPTION", description)
    mk(project, "UMBRELLA_PROJECT")

    if children:
        related = mk(project, "RELATED_PROJECTS")
        for child_accession in children:
            relation = mk(related, "RELATED_PROJECT")
            mk(relation, "CHILD_PROJECT", accession=child_accession)
    return project_set


def build_submission(action):
    root = ET.Element("SUBMISSION")
    actions = mk(root, "ACTIONS")
    action_node = mk(actions, "ACTION")
    mk(action_node, action.upper())
    return root


def build_xml(action, alias, title, description, children, accession=None):
    root = ET.Element("WEBIN_SUBMISSION")
    submission = build_submission(action)
    project_set = build_project(alias, title, description, children, accession)
    root.append(submission)
    root.append(project_set)
    return ET.ElementTree(root)

def write_tree(root, path):
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--action", choices=["add", "modify"], required=True)
    parser.add_argument("--alias", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--description", required=True)
    parser.add_argument("--project-accessions", nargs="+", default=[])
    parser.add_argument("--analysis-accessions", nargs="+", default=[])
    parser.add_argument("--umbrella-accession")
    parser.add_argument("--format", choices=["webin", "dropbox"], default="webin")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    children = child_accessions(args.project_accessions) if args.project_accessions else []
    if args.action == "modify" and not children:
        parser.error("--action modify requires project poll output containing a child PRJ accession")
    if args.action == "modify" and args.analysis_accessions:
        validate_analysis(args.analysis_accessions)

    os.makedirs(args.outdir, exist_ok=True)
    output = os.path.join(args.outdir, "umbrella_project.xml")
    if args.format == "dropbox":
        write_tree(
            build_project(args.alias, args.title, args.description, children, args.umbrella_accession),
            output,
        )
        write_tree(build_submission(args.action), os.path.join(args.outdir, "umbrella_submission.xml"))
    else:
        build_xml(args.action, args.alias, args.title, args.description, children, args.umbrella_accession).write(
            output, encoding="UTF-8", xml_declaration=True
        )
    print(output)


if __name__ == "__main__":
    main()
