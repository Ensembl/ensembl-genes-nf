#!/usr/bin/env python3
import argparse, os, sys, csv, datetime, xml.etree.ElementTree as ET

NS = {}

def mk_text(parent, tag, text=None, **attrib):
    el = ET.SubElement(parent, tag, attrib)
    if text:
        el.text = text
    return el


def build_analysis_xml(row, remote_path, md5, mode='READ_ALIGNMENT'):
    # Allowed analysis types per current ENA docs: READ_ALIGNMENT or REFERENCE_ALIGNMENT
    analysis_type = mode

    root = ET.Element('ANALYSIS_SET')
    analysis = mk_text(root, 'ANALYSIS', alias=row.get('analysis_alias') or f"auto-{os.path.basename(remote_path)}")

    mk_text(analysis, 'TITLE', row.get('title') or os.path.basename(remote_path))
    mk_text(analysis, 'DESCRIPTION', row.get('description') or 'Alignment of public runs')

    study = row.get('study')
    if study:
        mk_text(analysis, 'STUDY_REF', accession=study) if study.upper().startswith(('PRJ', 'ERP', 'SRP', 'DRP')) else mk_text(analysis, 'STUDY_REF', refname=study)

    sample = row.get('sample_accession')
    if sample:
        mk_text(analysis, 'SAMPLE_REF', accession=sample)

    # Link to runs
    runs = [r.strip() for r in (row.get('run_accessions') or '').split(',') if r.strip()]
    if runs:
        dset = mk_text(analysis, 'DATA_SET')
        for r in runs:
            mk_text(dset, 'RUN_REF', accession=r)

    files = mk_text(analysis, 'FILES')
    ftype = (row.get('file_type') or '').lower()
    filetype_attr = 'bam' if ftype == 'bam' else 'cram'
    mk_text(files, 'FILE', filename=remote_path, filetype=filetype_attr, checksum_method='MD5', checksum=md5)

    atype = mk_text(analysis, 'ANALYSIS_TYPE')
    assembly = row.get('assembly_accession')
    if analysis_type == 'READ_ALIGNMENT':
        ra = mk_text(atype, 'READ_ALIGNMENT')
        if assembly:
            # ENA expects ASSEMBLY/STANDARD accession
            asm = mk_text(ra, 'ASSEMBLY')
            mk_text(asm, 'STANDARD', accession=assembly)
    else:
        ra = mk_text(atype, 'REFERENCE_ALIGNMENT')
        if assembly:
            asm = mk_text(ra, 'ASSEMBLY')
            mk_text(asm, 'STANDARD', accession=assembly)
        # Optional sequences
        seqs = [s.strip() for s in (row.get('ref_seqs') or '').split(',') if s.strip()]
        for s in seqs:
            mk_text(ra, 'SEQUENCE', accession=s)

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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest-row', required=True, help='Path to a one-row TSV with headers')
    p.add_argument('--remote-path', required=True)
    p.add_argument('--md5', required=True)
    p.add_argument('--analysis-type', default='READ_ALIGNMENT', choices=['READ_ALIGNMENT','REFERENCE_ALIGNMENT'])
    p.add_argument('--hold-until', default=None)
    p.add_argument('--outdir', required=True)
    args = p.parse_args()

    with open(args.manifest_row, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        row = next(reader)

    md5val = None
    with open(args.md5) as m:
        md5val = m.read().strip()

    os.makedirs(args.outdir, exist_ok=True)
    analysis_tree = build_analysis_xml(row, args.remote_path, md5val, args.analysis_type)
    analysis_path = os.path.join(args.outdir, 'analysis.xml')
    analysis_tree.write(analysis_path, encoding='UTF-8', xml_declaration=True)

    submission_tree = build_submission_xml(args.hold_until)
    submission_path = os.path.join(args.outdir, 'submission.xml')
    submission_tree.write(submission_path, encoding='UTF-8', xml_declaration=True)

    print(analysis_path)
    print(submission_path)

if __name__ == '__main__':
    main()
