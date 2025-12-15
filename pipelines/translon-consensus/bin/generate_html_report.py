#!/usr/bin/env python3
"""
Generate comprehensive HTML report for translon consensus analysis
with detailed per-tool disagreement and performance summaries.

This script consumes *.consensus.tsv outputs produced by the ncORF
consensus analysis pipeline and aggregates them into a single
interactive HTML report.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime
import pandas as pd
from jinja2 import Template

# ---------------------------------------------------------------------
# Disagreement classification
# ---------------------------------------------------------------------

def classify_disagreement(start_matches: int, end_matches: int) -> str:
    """
    Classify disagreement type based on start/end matches.
    """
    if start_matches > 0 and end_matches > 0:
        return 'EXACT_MATCH'
    if start_matches > 0 and end_matches == 0:
        return 'START_ONLY'
    if start_matches == 0 and end_matches > 0:
        return 'END_ONLY'
    return 'BOTH_MISMATCH'


# ---------------------------------------------------------------------
# Template loading
# ---------------------------------------------------------------------

def load_template(template_path=None):
    if template_path:
        template_path = Path(template_path)
        if not template_path.exists():
            raise FileNotFoundError(f"Template not found: {template_path}")
        return Template(template_path.read_text())

    script_dir = Path(__file__).parent
    default = script_dir / '..' / 'assets' / 'report_template.html'
    if not default.exists():
        raise FileNotFoundError(
            f"Default template not found at {default}. "
            "Provide --template explicitly."
        )
    return Template(default.read_text())


# ---------------------------------------------------------------------
# Core parsing logic
# ---------------------------------------------------------------------

def parse_results_directory(results_dir: Path):
    """
    Parse all *.consensus.tsv files into a rich data structure suitable
    for reporting.
    """
    results_dir = Path(results_dir)
    samples = {}

    # Search for both *.consensus.tsv AND *.tsv that are actually consensus files
    all_tsv_files = list(results_dir.rglob('*.tsv'))
    
    for tsv_file in all_tsv_files:
        # Resolve symlinks to get actual filename
        actual_file = tsv_file.resolve()
        actual_name = actual_file.name
        
        # Skip if not a consensus file
        if '.consensus.tsv' not in actual_name:
            continue
            
        # Parse the actual filename: ncORFs_{SAMPLE}.{TOOL}.consensus.tsv
        name_without_ext = actual_name.replace('.consensus.tsv', '')
        
        # Remove 'ncORFs_' prefix if present
        if name_without_ext.startswith('ncORFs_'):
            name_without_ext = name_without_ext[7:]  # Remove 'ncORFs_'
        
        # Split on last dot to get sample and tool
        parts = name_without_ext.rsplit('.', 1)
        
        if len(parts) < 2:
            print(f"WARNING: Skipping {actual_name} - unexpected format", file=sys.stderr)
            continue

        sample_name, query_tool = parts[0], parts[1]

        if sample_name not in samples:
            samples[sample_name] = {
                'name': sample_name,
                'features': [],
                'tools': set(),
                'tool_disagreements': defaultdict(Counter),
                'pairwise': defaultdict(lambda: defaultdict(Counter)),
                'context_disagreements': defaultdict(lambda: defaultdict(Counter))
            }

        df = pd.read_csv(tsv_file, sep='\t')

        match_cols = [c for c in df.columns if c.endswith('_both_matches')]
        other_tools = [c.replace('_both_matches', '') for c in match_cols]

        samples[sample_name]['tools'].add(query_tool)
        samples[sample_name]['tools'].update(other_tools)

        for _, row in df.iterrows():
            feature = {
                'chr': row['chr'],
                'start': int(row['start_pos']),
                'end': int(row['end_pos']),
                'strand': row.get('strand', ''),
                'query_tool': query_tool,
                'feature_name': row.get('feature_name', 'N/A'),
                'rna_biotype': row.get('rna_biotype', 'intergenic'),
                'cds_context': row.get('cds_context', 'intergenic'),
                'gene_name': row.get('gene_name', 'intergenic'),
                'ucsc_url': row.get('ucsc_url', ''),
                'tool_matches': {}
            }

            for tool in other_tools:
                s = int(row.get(f'{tool}_start_matches', 0))
                e = int(row.get(f'{tool}_end_matches', 0))
                cls = classify_disagreement(s, e)

                feature['tool_matches'][tool] = {
                    'start': s,
                    'end': e,
                    'class': cls
                }

                # Aggregate per-tool disagreement
                samples[sample_name]['tool_disagreements'][tool][cls] += 1

                # Pairwise (directional)
                samples[sample_name]['pairwise'][query_tool][tool][cls] += 1

                # Context stratification
                ctx = feature['cds_context']
                samples[sample_name]['context_disagreements'][tool][ctx][cls] += 1

            samples[sample_name]['features'].append(feature)

    # Finalize tool lists
    for s in samples.values():
        s['tools'] = sorted(s['tools'])

    return {
        'samples': sorted(samples.values(), key=lambda x: x['name'])
    }


# ---------------------------------------------------------------------
# Signature summaries
# ---------------------------------------------------------------------

def compute_tool_signatures(sample):
    """
    Compute compact performance signatures per tool.
    """
    signatures = {}

    for tool, counts in sample['tool_disagreements'].items():
        total = sum(counts.values())
        if total == 0:
            continue

        exact = counts.get('EXACT_MATCH', 0)
        start_only = counts.get('START_ONLY', 0)
        end_only = counts.get('END_ONLY', 0)

        signatures[tool] = {
            'exact_match_rate': exact / total,
            'start_bias': (start_only - end_only) / total,
            'disagreement_rate': 1 - (exact / total),
            'raw_counts': dict(counts)
        }

    return signatures


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate HTML report from translon consensus results'
    )

    parser.add_argument('-i', '--results-dir', required=True)
    parser.add_argument('-o', '--output', default='translon_consensus_report.html')
    parser.add_argument('-n', '--run-name', default='Translon Consensus Analysis')
    parser.add_argument('-u', '--ucsc-session-url', default='https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38')
    parser.add_argument('-t', '--template', help='Custom HTML template')

    args = parser.parse_args()

    parsed = parse_results_directory(Path(args.results_dir))
    samples = parsed['samples']

    for s in samples:
        s['tool_signatures'] = compute_tool_signatures(s)

    total_features = sum(len(s['features']) for s in samples)

    context = {
        'run_name': args.run_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'ucsc_session_url': args.ucsc_session_url,
        'samples': samples,
        'total_samples': len(samples),
        'total_features': total_features
    }

    template = load_template(args.template)
    html = template.render(**context)

    output_path = Path(args.output)
    output_path.write_text(html)

    print(f"\n✅ HTML report generated: {output_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
