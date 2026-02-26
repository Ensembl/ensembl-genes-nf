#!/usr/bin/env python3
import argparse, os, sys, csv, re, datetime, xml.etree.ElementTree as ET

NS = {}

def mk_text(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, attrib)
    if text:
        el.text = text
    return el


def build_analysis_xml(row, remote_path, md5, mode='REFERENCE_ALIGNMENT'):
    root = ET.Element('ANALYSIS_SET')
    analysis = mk_text(root, 'ANALYSIS', alias=row.get('analysis_alias') or f"auto-{os.path.basename(remote_path)}")

    mk_text(analysis, 'TITLE', row.get('title') or os.path.basename(remote_path))
    mk_text(analysis, 'DESCRIPTION', row.get('description') or 'Alignment of public runs')

    study = row.get('study')
    if study:
        mk_text(analysis, 'STUDY_REF', accession=study) if study.upper().startswith(('PRJ', 'ERP', 'SRP', 'DRP')) else mk_text(analysis, 'STUDY_REF', refname=study)

    # One or more samples can be associated
    samples = [s.strip() for s in (row.get('sample_accession') or '').split(',') if s.strip()]
    for s in samples:
        mk_text(analysis, 'SAMPLE_REF', accession=s)

    # Link to runs/experiments directly under ANALYSIS (ENA schema)
    # Combine comma-separated list and optional file of run IDs
    runs = []
    inline_runs = [r.strip() for r in (row.get('run_accessions') or '').split(',') if r.strip()]
    runs.extend(inline_runs)

    run_list_path = (row.get('run_list_path') or '').strip()
    if run_list_path:
        try:
            with open(run_list_path) as rl:
                for line in rl:
                    # accept tokens separated by whitespace/commas/tabs; pick tokens that look like ERR/SRR/DRR
                    for tok in re.split(r'[\s,\t]+', line.strip()):
                        if re.match(r'^(ERR|SRR|DRR)\d{3,}$', tok, flags=re.IGNORECASE):
                            runs.append(tok)
        except FileNotFoundError:
            print(f"WARNING: run_list_path not found: {run_list_path}", file=sys.stderr)

    # de-duplicate while preserving order
    seen = set()
    runs = [x for x in runs if not (x in seen or seen.add(x))]

    if not row.get('__omit_run_refs__'):
        for r in runs:
            mk_text(analysis, 'RUN_REF', accession=r)

    exps = [e.strip() for e in (row.get('experiment_accessions') or '').split(',') if e.strip()]
    for e in exps:
        mk_text(analysis, 'EXPERIMENT_REF', accession=e)

    # ANALYSIS_TYPE followed by FILES
    atype = mk_text(analysis, 'ANALYSIS_TYPE')
    assembly = row.get('assembly_accession')
    analysis_type_tag = (mode or 'REFERENCE_ALIGNMENT').strip().upper()
    if analysis_type_tag == 'READ_ALIGNMENT':
        print('WARNING: READ_ALIGNMENT not supported by current ENA schema; mapping to REFERENCE_ALIGNMENT', file=sys.stderr)
        analysis_type_tag = 'REFERENCE_ALIGNMENT'
    if analysis_type_tag not in ('REFERENCE_ALIGNMENT',):
        print(f"WARNING: Unsupported analysis type '{mode}', defaulting to REFERENCE_ALIGNMENT", file=sys.stderr)
        analysis_type_tag = 'REFERENCE_ALIGNMENT'
    ra = mk_text(atype, analysis_type_tag)
    if assembly:
        asm = mk_text(ra, 'ASSEMBLY')
        mk_text(asm, 'STANDARD', accession=assembly)
    seqs = [s.strip() for s in (row.get('ref_seqs') or '').split(',') if s.strip()]
    for s in seqs:
        mk_text(ra, 'SEQUENCE', accession=s)

    files = mk_text(analysis, 'FILES')
    ftype = (row.get('file_type') or '').lower()
    filetype_attr = 'bam' if ftype == 'bam' else 'cram'
    mk_text(files, 'FILE', filename=remote_path, filetype=filetype_attr, checksum_method='MD5', checksum=md5)

    # Optional: ANALYSIS_LINKS (URL_LINK only for now) and ANALYSIS_ATTRIBUTES
    links_val = (row.get('analysis_links') or '').strip()
    if links_val:
        links_el = mk_text(analysis, 'ANALYSIS_LINKS')
        # Format: label|url; label2|url2
        for part in [p.strip() for p in links_val.split(';') if p.strip()]:
            try:
                label, url = [x.strip() for x in part.split('|', 1)]
            except ValueError:
                # If only URL provided, use URL as label
                label, url = part, part
            al = mk_text(links_el, 'ANALYSIS_LINK')
            ul = mk_text(al, 'URL_LINK')
            mk_text(ul, 'LABEL', label)
            mk_text(ul, 'URL', url)

    # Attributes from either a semicolon list or attr_* columns
    attrs = []
    attrs_list = (row.get('analysis_attributes') or '').strip()
    if attrs_list:
        # Format: key=value; key2=value2
        for part in [p.strip() for p in attrs_list.split(';') if p.strip()]:
            if '=' in part:
                k, v = [x.strip() for x in part.split('=', 1)]
                if k and v:
                    attrs.append((k, v))
    # Also accept any columns named attr_* as attributes
    for k, v in row.items():
        if k.startswith('attr_') and (v or '').strip():
            tag = k[len('attr_'):].strip()
            if tag:
                attrs.append((tag, v.strip()))

    if attrs:
        attrs_el = mk_text(analysis, 'ANALYSIS_ATTRIBUTES')
        for k, v in attrs:
            ael = mk_text(attrs_el, 'ANALYSIS_ATTRIBUTE')
            mk_text(ael, 'TAG', k)
            mk_text(ael, 'VALUE', v)

    return ET.ElementTree(root)


def build_submission_xml(hold_until=None):
    sub = ET.Element('SUBMISSION')
    acts = mk_text(sub, 'ACTIONS')
    if hold_until:
        mk_text(acts, 'ACTION')
        act = acts[-1]
        mk_text(act, 'HOLD', HoldUntilDate=hold_until)
    mk_text(acts, 'ACTION')
    act = acts[-1]
    mk_text(act, 'ADD')
    return ET.ElementTree(sub)


def build_webin_submission_xml(analysis_tree, submission_tree):
    root = ET.Element('WEBIN_SUBMISSION')
    root.append(submission_tree.getroot())
    root.append(analysis_tree.getroot())
    return ET.ElementTree(root)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest-row', required=True, help='Path to a one-row TSV with headers')
    p.add_argument('--remote-path', required=True)
    p.add_argument('--md5', required=True)
    p.add_argument('--analysis-type', default='REFERENCE_ALIGNMENT', choices=['READ_ALIGNMENT', 'REFERENCE_ALIGNMENT'])
    p.add_argument('--omit-run-refs', action='store_true', help='Omit RUN_REF elements (useful in ENA TEST)')
    p.add_argument('--hold-until', default=None)
    p.add_argument('--outdir', required=True)
    args = p.parse_args()

    with open(args.manifest_row, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        row = next(reader)
        if args.omit_run_refs:
            # pass an internal flag via the row dict to avoid changing function signature
            row['__omit_run_refs__'] = '1'

    with open(args.md5) as m:
        md5val = m.read().split()[0]

    os.makedirs(args.outdir, exist_ok=True)

    analysis_tree = build_analysis_xml(row, args.remote_path, md5val, args.analysis_type)
    analysis_path = os.path.join(args.outdir, 'analysis.xml')
    analysis_tree.write(analysis_path, encoding='UTF-8', xml_declaration=True)

    submission_tree = build_submission_xml(args.hold_until)
    submission_path = os.path.join(args.outdir, 'submission.xml')
    submission_tree.write(submission_path, encoding='UTF-8', xml_declaration=True)

    webin_tree = build_webin_submission_xml(analysis_tree, submission_tree)
    webin_path = os.path.join(args.outdir, 'webin_submission.xml')
    webin_tree.write(webin_path, encoding='UTF-8', xml_declaration=True)

    print(webin_path)

if __name__ == '__main__':
    main()
