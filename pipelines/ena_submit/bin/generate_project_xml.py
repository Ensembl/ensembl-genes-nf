#!/usr/bin/env python3
import argparse, os, sys, xml.etree.ElementTree as ET


def mk(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, attrib)
    if text is not None:
        el.text = text
    return el


def validate_text(name, value, minimum, maximum):
    length = len(value or '')
    if length < minimum or length > maximum:
        raise ValueError(f"{name} must contain {minimum}-{maximum} characters; got {length}")


def build_webin_project(alias, name, title, description, hold_until=None, scientific_name=None, taxon_id=None):
    validate_text('TITLE', title, 20, 250)
    validate_text('DESCRIPTION', description, 20, 4000)
    root = ET.Element('WEBIN')

    # SUBMISSION with ADD (and optional HOLD)
    submission_set = mk(root, 'SUBMISSION_SET')
    sub = mk(submission_set, 'SUBMISSION')
    acts = mk(sub, 'ACTIONS')
    mk(acts, 'ACTION'); mk(acts[-1], 'ADD')
    if hold_until:
        mk(acts, 'ACTION'); mk(acts[-1], 'HOLD', HoldUntilDate=hold_until)

    # PROJECT_SET
    prj_set = mk(root, 'PROJECT_SET')
    prj = mk(prj_set, 'PROJECT', alias=alias)
    if name:
        mk(prj, 'NAME', name)
    if title:
        mk(prj, 'TITLE', title)
    if description:
        mk(prj, 'DESCRIPTION', description)
    sub_prj = mk(prj, 'SUBMISSION_PROJECT')
    mk(sub_prj, 'SEQUENCING_PROJECT')
    if scientific_name or taxon_id:
        organism = mk(sub_prj, 'ORGANISM')
        if taxon_id:
            mk(organism, 'TAXON_ID', taxon_id)
        if scientific_name:
            mk(organism, 'SCIENTIFIC_NAME', scientific_name)

    return ET.ElementTree(root)


def main():
    ap = argparse.ArgumentParser(description='Generate Webin XML for registering a Project (Study)')
    ap.add_argument('--alias', required=True)
    ap.add_argument('--name', default='')
    ap.add_argument('--title', required=True)
    ap.add_argument('--description', required=True)
    ap.add_argument('--scientific-name', default='')
    ap.add_argument('--taxon-id', default='')
    ap.add_argument('--hold-until', default=None)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tree = build_webin_project(args.alias, args.name, args.title, args.description, args.hold_until, args.scientific_name, args.taxon_id)
    out_path = os.path.join(args.outdir, 'webin_project.xml')
    tree.write(out_path, encoding='UTF-8', xml_declaration=True)
    print(out_path)


if __name__ == '__main__':
    main()
