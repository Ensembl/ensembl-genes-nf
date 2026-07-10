#!/usr/bin/env python3
import argparse, csv, os, re, sys, xml.etree.ElementTree as ET

NS = {}

def mk_text(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, attrib)
    if text:
        el.text = text
    return el


def split_values(value):
    return [x.strip() for x in (value or '').split(',') if x.strip()]


def unique(values):
    seen = set()
    out = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def read_md5(path):
    with open(path) as handle:
        return handle.read().split()[0]


def read_files_manifest(path):
    rows = []
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            if not row.get('remote_path'):
                raise RuntimeError(f"files manifest row missing remote_path in {path}")
            if not row.get('file_type'):
                raise RuntimeError(f"files manifest row missing file_type for {row.get('remote_path')}")
            if row.get('md5'):
                md5 = row['md5'].strip()
            elif row.get('md5_path'):
                md5 = read_md5(row['md5_path'].strip())
            else:
                raise RuntimeError(f"files manifest row missing md5/md5_path for {row.get('remote_path')}")
            row['md5'] = md5
            rows.append(row)
    if not rows:
        raise RuntimeError(f"No files found in files manifest: {path}")
    return rows


def build_analysis_xml(row, files, mode='REFERENCE_ALIGNMENT'):
    root = ET.Element('ANALYSIS_SET')
    analysis = mk_text(root, 'ANALYSIS', alias=row.get('analysis_alias') or "rnaseq_alignment_evidence")

    mk_text(analysis, 'TITLE', row.get('title') or row.get('analysis_alias') or 'RNA-seq alignment evidence')
    mk_text(analysis, 'DESCRIPTION', row.get('description') or 'Alignment of public runs')

    study = (row.get('study') or '').strip()
    if study:
        # Treat as accession only if it matches an INSDC study/project accession pattern
        # e.g., PRJEB12345, PRJNA12345, ERP123456, SRP123456, DRP123456
        if re.match(r'^(PRJ[EDNA][A-Z]*\d+|ER[PS]\d+|SRP\d+|DRP\d+)$', study, flags=re.IGNORECASE):
            mk_text(analysis, 'STUDY_REF', accession=study)
        else:
            mk_text(analysis, 'STUDY_REF', refname=study)

    samples = split_values(row.get('sample_accession'))
    for file_row in files:
        samples.extend(split_values(file_row.get('sample_accession')))
    for s in unique(samples):
        mk_text(analysis, 'SAMPLE_REF', accession=s)

    # Link to runs/experiments directly under ANALYSIS (ENA schema)
    # Combine comma-separated list and optional file of run IDs
    runs = []
    runs.extend(split_values(row.get('run_accessions')))
    for file_row in files:
        runs.extend(split_values(file_row.get('run_accession') or file_row.get('run_accessions')))

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

    if not row.get('__omit_run_refs__'):
        for r in unique(runs):
            mk_text(analysis, 'RUN_REF', accession=r)

    exps = split_values(row.get('experiment_accessions'))
    for file_row in files:
        exps.extend(split_values(file_row.get('experiment_accession') or file_row.get('experiment_accessions')))
    for e in unique(exps):
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

    # FILES must appear before ANALYSIS_LINKS / ANALYSIS_ATTRIBUTES
    files_el = mk_text(analysis, 'FILES')
    for file_row in files:
        ftype = (file_row.get('file_type') or '').lower()
        if ftype not in ('bam', 'cram'):
            raise RuntimeError(f"Unsupported file_type '{ftype}' for {file_row.get('remote_path')}")
        mk_text(
            files_el,
            'FILE',
            filename=file_row['remote_path'],
            filetype=ftype,
            checksum_method='MD5',
            checksum=file_row['md5'],
        )

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
    p.add_argument('--files-manifest', required=True, help='TSV with remote_path, file_type, md5 or md5_path')
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

    files = read_files_manifest(args.files_manifest)

    os.makedirs(args.outdir, exist_ok=True)

    analysis_tree = build_analysis_xml(row, files, args.analysis_type)
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
