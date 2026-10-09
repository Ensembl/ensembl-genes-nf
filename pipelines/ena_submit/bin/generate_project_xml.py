#!/usr/bin/env python3
import argparse, os, sys, xml.etree.ElementTree as ET


def mk(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, attrib)
    if text is not None:
        el.text = text
    return el


def build_webin_project(alias, name, title, description, hold_until=None):
    # Use WEBIN_SUBMISSION wrapper compatible with Webin v2 queue endpoint
    root = ET.Element('WEBIN_SUBMISSION')

    # SUBMISSION with ADD (and optional HOLD)
    sub = mk(root, 'SUBMISSION')
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

    return ET.ElementTree(root)


def main():
    ap = argparse.ArgumentParser(description='Generate Webin XML for registering a Project (Study)')
    ap.add_argument('--alias', required=True)
    ap.add_argument('--name', default='')
    ap.add_argument('--title', required=True)
    ap.add_argument('--description', required=True)
    ap.add_argument('--hold-until', default=None)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tree = build_webin_project(args.alias, args.name, args.title, args.description, args.hold_until)
    out_path = os.path.join(args.outdir, 'webin_project.xml')
    tree.write(out_path, encoding='UTF-8', xml_declaration=True)
    print(out_path)


if __name__ == '__main__':
    main()
