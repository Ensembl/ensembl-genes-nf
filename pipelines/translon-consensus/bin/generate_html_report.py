#!/usr/bin/env python3
"""
Generate comprehensive HTML report for translon consensus analysis.
Aggregates results from all samples into a single interactive HTML report.
"""

import argparse
import sys
import json
from pathlib import Path
from collections import defaultdict
import pandas as pd
from jinja2 import Template


def load_template(template_path=None):
    """
    Load the HTML template from assets directory or use default path.

    Args:
        template_path: Optional path to template file

    Returns:
        Jinja2 Template object
    """
    if template_path and Path(template_path).exists():
        template_file = Path(template_path)
    else:
        # Try to find template relative to this script
        script_dir = Path(__file__).parent
        template_file = script_dir.parent / 'assets' / 'report_template.html'

        if not template_file.exists():
            raise FileNotFoundError(
                f"Could not find template at {template_file}. "
                f"Please ensure report_template.html exists in the assets directory."
            )

    return Template(template_file.read_text())


def parse_results_directory(results_dir):
    """
    Parse all consensus results into a feature-centric format.

    Returns:
        dict with keys: samples (list of sample data with features)
    """
    results_dir = Path(results_dir)
    samples = {}

    # Find all consensus detail files
    for consensus_file in results_dir.rglob('*.consensus.tsv'):
        parts = consensus_file.stem.split('.')
        if len(parts) >= 3:  # sample.tool.consensus
            sample_name = parts[0]
            query_tool = parts[1]

            if sample_name not in samples:
                samples[sample_name] = {
                    'name': sample_name,
                    'tools': set(),
                    'features': [],
                    'consensus_stats': {}  # Will aggregate across tools
                }

            try:
                df = pd.read_csv(consensus_file, sep='\t')

                # Get other tool names from column names (those ending in _both_matches)
                match_cols = [col for col in df.columns if col.endswith('_both_matches')]
                other_tools = [col.replace('_both_matches', '') for col in match_cols]

                samples[sample_name]['tools'].add(query_tool)
                samples[sample_name]['tools'].update(other_tools)

                # Calculate consensus scores for all rows first
                df['consensus_score'] = df[[col for col in match_cols]].apply(
                    lambda row: sum(1 for val in row if val > 0), axis=1
                )

                # Group by consensus score and sample up to 500 from each level
                sampled_features = []
                consensus_stats = {}

                for score in sorted(df['consensus_score'].unique(), reverse=True):
                    score_group = df[df['consensus_score'] == score]
                    total_count = len(score_group)
                    sample_count = min(500, total_count)
                    sampled = score_group.head(sample_count)

                    sampled_features.append(sampled)
                    consensus_stats[int(score)] = {
                        'total': total_count,
                        'sampled': sample_count
                    }

                df_sorted = pd.concat(sampled_features, ignore_index=True)

                print(f"  {query_tool}: Sampled {len(df_sorted)} of {len(df)} features across {len(consensus_stats)} consensus levels", file=sys.stderr)

                # Aggregate consensus stats for this sample
                for score, stats in consensus_stats.items():
                    if score not in samples[sample_name]['consensus_stats']:
                        samples[sample_name]['consensus_stats'][score] = {
                            'total': 0,
                            'sampled': 0
                        }
                    samples[sample_name]['consensus_stats'][score]['total'] += stats['total']
                    samples[sample_name]['consensus_stats'][score]['sampled'] += stats['sampled']

                # Process features
                for _, row in df_sorted.iterrows():
                    consensus_count = row['consensus_score']
                    total_tools = len(match_cols) + 1  # +1 for query tool itself

                    feature = {
                        'sample': sample_name,
                        'query_tool': query_tool,
                        'chr': row['chr'],
                        'start': row['start_pos'],
                        'end': row['end_pos'],
                        'strand': row['strand'],
                        'name': row.get('feature_name', 'N/A'),
                        'consensus_score': int(consensus_count),
                        'total_tools': total_tools,
                        'ucsc_url': row.get('ucsc_url', ''),
                        'rna_biotype': row.get('rna_biotype', 'N/A'),
                        'cds_context': row.get('cds_context', 'N/A'),
                        'gene_name': row.get('gene_name', 'N/A'),
                        'tool_matches': {}
                    }

                    # Add per-tool match info
                    for tool in other_tools:
                        feature['tool_matches'][tool] = {
                            'start': int(row.get(f'{tool}_start_matches', 0)),
                            'end': int(row.get(f'{tool}_end_matches', 0)),
                            'both': int(row.get(f'{tool}_both_matches', 0))
                        }

                    samples[sample_name]['features'].append(feature)

            except Exception as e:
                print(f"Warning: Failed to parse {consensus_file}: {e}", file=sys.stderr)

    # Convert sets to sorted lists
    for sample_data in samples.values():
        sample_data['tools'] = sorted(list(sample_data['tools']))
        sample_data['tool_count'] = len(sample_data['tools'])

    return {
        'samples': sorted(samples.values(), key=lambda x: x['name'])
    }


def main():
    parser = argparse.ArgumentParser(
        description='Generate HTML report from translon consensus results'
    )

    parser.add_argument(
        '-i', '--results-dir',
        required=True,
        help='Directory containing consensus analysis results'
    )

    parser.add_argument(
        '-o', '--output',
        default='translon_consensus_report.html',
        help='Output HTML file (default: translon_consensus_report.html)'
    )

    parser.add_argument(
        '-n', '--run-name',
        default='Translon Consensus Analysis',
        help='Name of the pipeline run'
    )

    parser.add_argument(
        '-u', '--ucsc-session-url',
        default='https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38',
        help='UCSC Genome Browser session URL'
    )

    parser.add_argument(
        '-t', '--template',
        help='Path to custom HTML template file (default: uses assets/report_template.html)'
    )

    args = parser.parse_args()

    print(f"Parsing results from: {args.results_dir}", file=sys.stderr)
    parsed_data = parse_results_directory(args.results_dir)

    samples = parsed_data['samples']
    print(f"Found {len(samples)} samples", file=sys.stderr)

    total_features = sum(len(s['features']) for s in samples)
    print(f"Found {total_features} total features", file=sys.stderr)

    # Prepare template context
    from datetime import datetime
    context = {
        'run_name': args.run_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'ucsc_session_url': args.ucsc_session_url,
        'samples': samples,
        'total_samples': len(samples),
        'total_features': total_features
    }

    # Load and render template
    print(f"Loading HTML template...", file=sys.stderr)
    template = load_template(args.template)
    html_output = template.render(**context)

    # Write output
    output_path = Path(args.output)
    output_path.write_text(html_output)

    print(f"\n✅ HTML report generated: {output_path}", file=sys.stderr)
    print(f"   Open in browser: file://{output_path.absolute()}", file=sys.stderr)


if __name__ == '__main__':
    main()
