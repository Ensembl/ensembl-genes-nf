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


def normalise_attribute_key(key):
    return key[5:] if key.startswith('attr_') else key


ENA_SAMPLE_RE = re.compile(r'^(ERS|SRS|DRS)\d+$', re.IGNORECASE)
BIO_SAMPLE_RE = re.compile(r'^SAM(?:EA|D|N)\d+$', re.IGNORECASE)


def looks_like_sample_accession(value):
    """Return whether value is an ENA Sample accession usable in SAMPLE_REF."""
    return bool(ENA_SAMPLE_RE.match((value or '').strip()))


def read_md5(path):
    with open(path) as handle:
        return handle.read().split()[0]


def validate_md5(value, remote_path):
    if not re.fullmatch(r"[0-9a-fA-F]{32}", value or ""):
        raise RuntimeError(f"invalid MD5 checksum for {remote_path}: {value!r}")
    return value.lower()


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
            row['md5'] = validate_md5(md5, row['remote_path'])
            rows.append(row)
    if not rows:
        raise RuntimeError(f"No files found in files manifest: {path}")
    if len(rows) != 1:
        raise RuntimeError(
            f"Each ENA REFERENCE_ALIGNMENT analysis must contain exactly one BAM or CRAM; "
            f"found {len(rows)} in {path}"
        )
    return rows


def build_analysis_xml(row, files, mode='REFERENCE_ALIGNMENT'):
    root = ET.Element('ANALYSIS_SET')
    analysis = mk_text(root, 'ANALYSIS', alias=row.get('analysis_alias') or "rnaseq_alignment_evidence")

    mk_text(analysis, 'TITLE', row.get('title') or row.get('analysis_alias') or 'RNA-seq alignment evidence')
    mk_text(analysis, 'DESCRIPTION', row.get('description') or 'Alignment of public runs')

    study = (row.get('study') or '').strip()
    if not study:
        raise RuntimeError('REFERENCE_ALIGNMENT analysis requires a study')
    if study:
        # Treat as accession only if it matches an INSDC study/project accession pattern
        # e.g., PRJEB12345, PRJNA12345, ERP123456, SRP123456, DRP123456
        if re.match(r'^(PRJ[A-Z]+\d+|ERP\d+|SRP\d+|DRP\d+)$', study, flags=re.IGNORECASE):
            mk_text(analysis, 'STUDY_REF', accession=study)
        else:
            mk_text(analysis, 'STUDY_REF', refname=study)

    samples = split_values(row.get('sample_accession'))
    samples.extend(split_values(row.get('sample_accessions')))
    for file_row in files:
        samples.extend(split_values(file_row.get('sample_accession')))
        samples.extend(split_values(file_row.get('sample_accessions')))
    invalid_sample_refs = [s for s in unique(samples) if BIO_SAMPLE_RE.match(s)]
    if invalid_sample_refs:
        raise RuntimeError(
            "BioSample accession(s) cannot be emitted as ANALYSIS/SAMPLE_REF; "
            "resolve each source run to its ENA Sample accession (ERS/SRS/DRS): "
            + ", ".join(invalid_sample_refs)
        )
    sample_accessions = unique(samples)
    biosample_accessions = split_values(row.get('biosample_accessions'))
    for file_row in files:
        biosample_accessions.extend(split_values(file_row.get('biosample_accessions')))
    biosample_accessions = unique(biosample_accessions)
    for index, s in enumerate(sample_accessions):
        if looks_like_sample_accession(s):
            sample_ref = mk_text(analysis, 'SAMPLE_REF', accession=s)
            if len(sample_accessions) == len(biosample_accessions):
                identifiers = mk_text(sample_ref, 'IDENTIFIERS')
                mk_text(identifiers, 'PRIMARY_ID', s)
                mk_text(identifiers, 'EXTERNAL_ID', biosample_accessions[index], namespace='BioSample')
        elif s:
            print(f"WARNING: ignoring non-accession sample value '{s}'", file=sys.stderr)

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

    exps = split_values(row.get('experiment_accessions'))
    for file_row in files:
        exps.extend(split_values(file_row.get('experiment_accession') or file_row.get('experiment_accessions')))
    # SRA.analysis.xsd requires EXPERIMENT_REF before RUN_REF.
    for e in unique(exps):
        mk_text(analysis, 'EXPERIMENT_REF', accession=e)
    if not row.get('__omit_run_refs__'):
        for r in unique(runs):
            mk_text(analysis, 'RUN_REF', accession=r)

    # ANALYSIS_TYPE followed by FILES
    atype = mk_text(analysis, 'ANALYSIS_TYPE')
    assembly = (row.get('assembly_accession') or '').strip()
    ref_seqs = [s.strip() for s in (row.get('ref_seqs') or '').split(',') if s.strip()]
    if not assembly and not ref_seqs:
        raise RuntimeError(
            'REFERENCE_ALIGNMENT analysis requires assembly_accession or ref_seqs'
        )
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
    for s in ref_seqs:
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
                    attrs.append((normalise_attribute_key(k), v))
    # Also accept any columns named attr_* as attributes
    for k, v in row.items():
        if k.startswith('attr_') and (v or '').strip():
            tag = normalise_attribute_key(k).strip()
            if tag:
                attrs.append((tag, v.strip()))

    # These values are resolved from ENA or declared by the manifest and are
    # useful as searchable analysis-level metadata, while the authoritative
    # biological taxonomy remains on the referenced SAMPLE records.
    for key in ('scientific_name', 'species', 'taxon_id', 'aggregation_scope', 'metadata_source',
                'platform', 'source_fastq_count', 'alignment_software', 'alignment_software_version',
                'bam_sort_order', 'bam_pg_programs', 'source_metadata_count',
                'library_strategies', 'library_sources', 'instrument_platforms',
                'instrument_models', 'source_study_accessions',
                'secondary_study_accessions', 'secondary_sample_accessions',
                'biosample_accessions'):
        value = (row.get(key) or '').strip()
        if value and not any(existing == key for existing, _ in attrs):
            attrs.append((key, value))

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
    mk_text(acts, 'ACTION')
    act = acts[-1]
    mk_text(act, 'ADD')
    if hold_until:
        mk_text(acts, 'ACTION')
        act = acts[-1]
        mk_text(act, 'HOLD', HoldUntilDate=hold_until)
    return ET.ElementTree(sub)


def build_webin_submission_xml(analysis_tree, submission_tree):
    root = ET.Element('WEBIN')
    submission_set = mk_text(root, 'SUBMISSION_SET')
    submission_set.append(submission_tree.getroot())
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
