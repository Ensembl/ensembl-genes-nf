#!/usr/bin/env python3
import csv, argparse, os, sys, re, pathlib, json
from collections import defaultdict, Counter


def slug(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip())
    return re.sub(r"_+", "_", s).strip("_")


def load_map(path, key='taxon_id', val='species'):
    m = {}
    if not path:
        return m
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter='\t' if path.endswith('.tsv') else ',')
        for r in rd:
            m[str(r[key]).strip()] = r[val].strip()
    return m


def choose_group_key(row, group_by):
    parts = []
    for g in group_by:
        if g == 'assembly':
            parts.append(row.get('assembly_accession') or row.get('assembly') or '')
        elif g == 'taxon':
            parts.append(str(row.get('taxon_id') or ''))
        elif g == 'tissue':
            parts.append(row.get('tissue_prediction') or '')
        else:
            parts.append(row.get(g, ''))
    return tuple(parts)


def build_alias(species, assembly, release, group_label=None):
    base = f"{slug(species)}_{slug(assembly)}_{slug(release)}_aln"
    if group_label:
        base = f"{base}_{slug(group_label)}"
    return base


def main():
    ap = argparse.ArgumentParser(description='Build ENA submission manifest from a runs/tissue CSV')
    ap.add_argument('--input', required=True, help='Input CSV with columns including run_accession, sample_accession, assembly_accession, taxon_id, tissue_prediction')
    ap.add_argument('--study', required=True, help='ENA Study accession or alias (PRJ...) to assign to all rows, or use --study-map')
    ap.add_argument('--study-map', help='TSV/CSV mapping of assembly_accession or taxon_id to study; columns: key,study. Overrides --study per matching row.')
    ap.add_argument('--species-map', help='TSV/CSV mapping of taxon_id to species; columns: taxon_id,species')
    ap.add_argument('--release', default='Ensembl_release', help='Release label used in titles/descriptions (e.g., Ensembl_110)')
    ap.add_argument('--mode', choices=['per-run','merged'], default='per-run', help='per-run: one analysis per run; merged: one per group')
    ap.add_argument('--group-by', default='assembly', help='When --mode merged, comma list of grouping keys: assembly[,tissue|taxon]')
    ap.add_argument('--file-dir', required=True, help='Directory containing final merged BAM/CRAM files')
    ap.add_argument('--file-ext', default='cram', choices=['bam','cram'], help='File extension and file_type')
    ap.add_argument('--outdir', required=True, help='Output directory for manifest and run lists')
    ap.add_argument('--projects-out', help='Optional path to write projects.tsv derived from input')
    ap.add_argument('--links-prefix', help='Optional URL prefix where run lists will be hosted; used to populate analysis_links')
    ap.add_argument('--min-pass', action='store_true', help='Filter rows to those with passed_both==True when column exists')
    ap.add_argument('--priority-max', type=int, default=None, help='Filter rows to priority <= N when column exists')

    args = ap.parse_args()
    group_by = [g.strip() for g in args.group_by.split(',') if g.strip()]
    species_map = load_map(args.species_map, key='taxon_id', val='species')

    # study_map supports keys by assembly_accession or taxon_id in one file
    study_map = {}
    if args.study_map:
        with open(args.study_map) as fh:
            rd = csv.DictReader(fh, delimiter='\t' if args.study_map.endswith('.tsv') else ',')
            for r in rd:
                k = (r.get('assembly_accession') or r.get('taxon_id') or r.get('key') or '').strip()
                if k:
                    study_map[k] = r['study'].strip()

    os.makedirs(args.outdir, exist_ok=True)
    runs_dir = os.path.join(args.outdir, 'runs')
    os.makedirs(runs_dir, exist_ok=True)

    groups = defaultdict(list)
    projects = {}
    with open(args.input) as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            if args.min_pass and 'passed_both' in row and str(row['passed_both']).lower() not in ('true','1','yes'):
                continue
            if args.priority_max is not None and 'priority' in row:
                try:
                    if int(row['priority']) > args.priority_max:
                        continue
                except Exception:
                    pass
            if args.mode == 'merged':
                key = choose_group_key(row, group_by)
                groups[key].append(row)
            else:
                # per-run: use a synthetic key per row
                key = ('__per_run__', row.get('run_accession') or row.get('run') or '')
                groups[key].append(row)
            # derive per-assembly project alias (workflow also derives; this keeps manifest readable)
            assembly = row.get('assembly_accession') or row.get('assembly') or ''
            taxon_id = str(row.get('taxon_id') or '').strip()
            palias = f"prj_{slug(assembly)}" if assembly else (f"prj_taxon_{slug(taxon_id)}" if taxon_id else None)
            if palias and palias not in projects:
                projects[palias] = {
                    'alias': palias,
                    'name': palias,
                    'title': f"Annotation evidence project for {assembly or taxon_id}",
                    'description': "Per-assembly project auto-generated from runs CSV"
                }

    manifest_path = os.path.join(args.outdir, 'manifest.tsv')
    files_dir = os.path.join(args.outdir, 'files')
    os.makedirs(files_dir, exist_ok=True)
    manifest_cols = [
        'files_tsv', 'study', 'project_alias', 'umbrella_study', 'analysis_alias',
        'title', 'description', 'assembly_accession', 'last_geneset_update',
        'partial_release_label', 'species', 'taxon_id', 'analysis_links',
        'analysis_attributes', 'analysis_type', 'omit_run_refs_in_test',
    ]
    file_cols = [
        'file_path', 'file_type', 'remote_name', 'run_accession', 'run_accessions',
        'sample_accession', 'sample_accessions', 'experiment_accession',
    ]
    with open(manifest_path, 'w', newline='') as mf:
        manifest_writer = csv.DictWriter(mf, fieldnames=manifest_cols, delimiter='\t')
        manifest_writer.writeheader()

        for key, rows in groups.items():
            # Common derived values
            assembly = next((r.get('assembly_accession') for r in rows if r.get('assembly_accession')), '')
            taxon_id = str(next((r.get('taxon_id') for r in rows if r.get('taxon_id')), '')).strip()
            species = species_map.get(taxon_id) or (f"taxon_{taxon_id}" if taxon_id else 'unknown_species')

            group_outputs = []
            if args.mode == 'per-run':
                for r in rows:
                    ra = (r.get('run_accession') or r.get('run') or '').strip()
                    if not ra:
                        continue
                    sa = (r.get('sample_accession') or '').strip()
                    tissue = (r.get('tissue_prediction') or '').strip()
                    alias = slug(f"{ra}_{assembly}_{args.release}_aln")
                    remote_name = f"{alias}.{args.file_ext}"
                    file_path = os.path.join(args.file_dir, remote_name)
                    title = f"Run {ra} aligned to {assembly} for {species}, {args.release}"
                    desc = f"Alignment of run {ra}{' ('+tissue+')' if tissue else ''} to {assembly}; evidence for genome annotation {args.release}."
                    links = ''
                    attrs = [("reference_accession", assembly), ("source_runs_count", "1")]
                    if tissue:
                        attrs.append(("tissue_summary", tissue))
                        attrs.append(("tissue_method", "prediction from metadata"))
                    analysis_attributes = '; '.join([f"attr_{k}={v}" if not k.startswith('attr_') else f"{k}={v}" for k,v in attrs])
                    study = args.study
                    if study_map:
                        study = study_map.get(assembly) or study_map.get(taxon_id) or study
                    group_outputs.append({
                        'file_path': file_path,
                        'file_type': args.file_ext,
                        'remote_name': remote_name,
                        'run_accession': ra,
                        'run_accessions': ra,
                        'sample_accession': sa,
                        'sample_accessions': sa,
                        'experiment_accession': (r.get('experiment_accession') or '').strip(),
                    })
                    group_manifest = {
                        'analysis_alias': alias, 'title': title, 'description': desc,
                        'analysis_links': links, 'analysis_attributes': analysis_attributes,
                    }
                    group_outputs[-1]['_manifest'] = group_manifest
            else:
                tissue = None
                if 'tissue' in group_by or 'tissue_prediction' in group_by:
                    tissue = key[group_by.index('tissue') if 'tissue' in group_by else group_by.index('tissue_prediction')]
                group_label = tissue if tissue else None
                alias = build_alias(species, assembly, args.release, group_label)
                # run list and counts
                runs = []
                sample_counts = Counter()
                tissue_counts = Counter()
                for r in rows:
                    ra = (r.get('run_accession') or r.get('run') or '').strip()
                    if ra:
                        runs.append(ra)
                    sa = (r.get('sample_accession') or '').strip()
                    if sa:
                        sample_counts[sa] += 1
                    tp = (r.get('tissue_prediction') or '').strip()
                    if tp:
                        tissue_counts[tp] += 1
                seen = set()
                run_list = [x for x in runs if not (x in seen or seen.add(x))]
                run_list_rel = f"runs/{alias}.tsv"
                run_list_path = os.path.join(args.outdir, run_list_rel)
                with open(run_list_path, 'w') as rl:
                    rl.write('run\tsample\ttissue\tstudy\n')
                    for r in rows:
                        ra = (r.get('run_accession') or r.get('run') or '').strip()
                        if not ra:
                            continue
                        sa = (r.get('sample_accession') or '').strip()
                        tp = (r.get('tissue_prediction') or '').strip()
                        st = (r.get('study') or r.get('study_accession') or '').strip()
                        rl.write(f"{ra}\t{sa}\t{tp}\t{st}\n")
                # links/attrs
                links = f"Runs TSV|{args.links_prefix.rstrip('/')}/{alias}.tsv" if args.links_prefix else ''
                attrs = [("reference_accession", assembly), ("source_runs_count", str(len(run_list)))]
                if tissue_counts:
                    parts = [f"{k}({v})" for k, v in tissue_counts.most_common(12)]
                    attrs.append(("tissue_summary", ", ".join(parts)))
                    attrs.append(("tissue_method", "prediction from metadata"))
                analysis_attributes = '; '.join([f"attr_{k}={v}" if not k.startswith('attr_') else f"{k}={v}" for k,v in attrs])
                study = args.study
                if study_map:
                    study = study_map.get(assembly) or study_map.get(taxon_id) or study
                title = f"Merged {'%s ' % tissue if tissue else ''}reads aligned to {assembly} for {species}, {args.release}"
                desc = f"Merged {len(run_list)} public runs across multiple BioSamples; evidence for genome annotation {args.release}."
                remote_name = f"{alias}.{args.file_ext}"
                file_path = os.path.join(args.file_dir, remote_name)
                group_outputs = [{
                    'file_path': file_path,
                    'file_type': args.file_ext,
                    'remote_name': remote_name,
                    'run_accessions': ','.join(run_list),
                    'sample_accessions': ','.join(sorted(sample_counts)),
                    '_manifest': {'analysis_alias': alias, 'title': title,
                                  'description': desc, 'analysis_links': links,
                                  'analysis_attributes': analysis_attributes},
                }]

            if not group_outputs:
                continue
            manifest_info = group_outputs[0].pop('_manifest')
            files_tsv = os.path.join(files_dir, f"{manifest_info['analysis_alias']}.tsv")
            with open(files_tsv, 'w', newline='') as ff:
                file_writer = csv.DictWriter(ff, fieldnames=file_cols, delimiter='\t')
                file_writer.writeheader()
                file_writer.writerows(group_outputs)
            project_alias = f"prj_{slug(assembly)}" if assembly else f"prj_taxon_{slug(taxon_id)}"
            manifest_writer.writerow({
                'files_tsv': files_tsv, 'study': study, 'project_alias': project_alias,
                'analysis_alias': manifest_info['analysis_alias'],
                'title': manifest_info['title'], 'description': manifest_info['description'],
                'assembly_accession': assembly, 'last_geneset_update': args.release,
                'partial_release_label': args.release, 'species': species, 'taxon_id': taxon_id,
                'analysis_links': manifest_info['analysis_links'],
                'analysis_attributes': manifest_info['analysis_attributes'],
                'analysis_type': 'REFERENCE_ALIGNMENT', 'omit_run_refs_in_test': 'true',
            })

    print(manifest_path)

    if args.projects_out and projects:
        with open(args.projects_out, 'w', newline='') as pf:
            cols = ['alias','name','title','description']
            w = csv.DictWriter(pf, fieldnames=cols, delimiter='\t')
            w.writeheader()
            for p in projects.values():
                w.writerow(p)


if __name__ == '__main__':
    main()
