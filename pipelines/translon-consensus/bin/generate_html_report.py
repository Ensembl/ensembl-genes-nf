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
    Parse all consensus results from a directory.

    Returns:
        dict with keys: samples_with_consensus, skipped_samples, all_tools
    """
    results_dir = Path(results_dir)

    samples_data = defaultdict(lambda: {
        'name': '',
        'tool_count': 0,
        'tools': set(),
        'summary_file': None,
        'summary_data': None,
        'consensus_files': {},
        'consensus_data': {}
    })

    # Find all summary files
    for summary_file in results_dir.rglob('*.summary.tsv'):
        sample_name = summary_file.stem.replace('.summary', '')

        try:
            summary_df = pd.read_csv(summary_file, sep='\t')
            samples_data[sample_name]['summary_file'] = str(summary_file)
            samples_data[sample_name]['summary_data'] = summary_df.to_dict('records')
            samples_data[sample_name]['name'] = sample_name
            samples_data[sample_name]['tools'] = set(summary_df['tool'].tolist())
            samples_data[sample_name]['tool_count'] = len(samples_data[sample_name]['tools'])
        except Exception as e:
            print(f"Warning: Failed to parse {summary_file}: {e}", file=sys.stderr)

    # Find all consensus detail files
    for consensus_file in results_dir.rglob('*.consensus.tsv'):
        parts = consensus_file.stem.split('.')
        if len(parts) >= 3:  # sample.tool.consensus
            sample_name = parts[0]
            tool = parts[1]

            try:
                consensus_df = pd.read_csv(consensus_file, sep='\t')
                # Limit to first 100 rows for HTML (performance)
                consensus_df_limited = consensus_df.head(100)

                samples_data[sample_name]['consensus_files'][tool] = str(consensus_file)
                samples_data[sample_name]['consensus_data'][tool] = consensus_df_limited.to_dict('records')
            except Exception as e:
                print(f"Warning: Failed to parse {consensus_file}: {e}", file=sys.stderr)

    # Separate multi-tool vs single-tool samples
    samples_with_consensus = []
    skipped_samples = []
    all_tools = set()

    for sample_name, data in samples_data.items():
        data['tools'] = sorted(list(data['tools']))
        all_tools.update(data['tools'])

        if data['tool_count'] >= 2:
            samples_with_consensus.append(data)
        else:
            skipped_samples.append(data)

    return {
        'samples_with_consensus': sorted(samples_with_consensus, key=lambda x: x['name']),
        'skipped_samples': sorted(skipped_samples, key=lambda x: x['name']),
        'all_tools': sorted(list(all_tools))
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

    print(f"Found {len(parsed_data['samples_with_consensus'])} samples with consensus", file=sys.stderr)
    print(f"Found {len(parsed_data['skipped_samples'])} single-tool samples", file=sys.stderr)

    # Calculate tool usage counts
    tool_counts_dict = defaultdict(int)
    for sample in parsed_data['samples_with_consensus'] + parsed_data['skipped_samples']:
        for tool in sample['tools']:
            tool_counts_dict[tool] += 1

    tool_names = sorted(tool_counts_dict.keys())
    tool_counts = [tool_counts_dict[t] for t in tool_names]

    # Prepare template context
    from datetime import datetime
    context = {
        'run_name': args.run_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'ucsc_session_url': args.ucsc_session_url,
        'total_samples': len(parsed_data['samples_with_consensus']) + len(parsed_data['skipped_samples']),
        'multi_tool_samples': len(parsed_data['samples_with_consensus']),
        'single_tool_samples': len(parsed_data['skipped_samples']),
        'unique_tools': len(parsed_data['all_tools']),
        'all_tools': parsed_data['all_tools'],
        'tool_names': tool_names,
        'tool_counts': tool_counts,
        'samples_with_consensus': parsed_data['samples_with_consensus'],
        'skipped_samples': parsed_data['skipped_samples']
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
