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
import numpy as np
from jinja2 import Template

# ---------------------------------------------------------------------
# Type conversion utilities
# ---------------------------------------------------------------------

def convert_numpy_types(obj):
    """
    Recursively convert numpy types to Python native types for JSON serialization.
    """
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy_types(item) for item in obj]
    elif isinstance(obj, Counter):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, defaultdict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    return obj

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

def parse_results_directory(results_dir: Path, max_examples_per_category: int = 100):
    """
    Parse all *.consensus.tsv files into a rich data structure suitable
    for reporting. Uses stratified sampling to limit examples per disagreement category.

    Args:
        results_dir: Directory containing consensus TSV files
        max_examples_per_category: Maximum number of feature examples to keep per disagreement type
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
                'features_by_consensus': defaultdict(list),  # Stratify by consensus score
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

        # OPTIMIZATION: Calculate consensus scores vectorized
        consensus_scores = np.zeros(len(df), dtype=int)
        for tool in other_tools:
            start_col = f'{tool}_start_matches'
            end_col = f'{tool}_end_matches'
            # EXACT_MATCH = both start and end match
            exact_match = (df[start_col] > 0) & (df[end_col] > 0)
            consensus_scores += exact_match.astype(int)

        # Add consensus_score column
        df['consensus_score'] = consensus_scores

        # OPTIMIZATION: Aggregate disagreements vectorized before sampling
        for tool in other_tools:
            start_col = f'{tool}_start_matches'
            end_col = f'{tool}_end_matches'

            # Vectorized classification
            has_start = df[start_col] > 0
            has_end = df[end_col] > 0

            exact = (has_start & has_end).sum()
            start_only = (has_start & ~has_end).sum()
            end_only = (~has_start & has_end).sum()
            neither = (~has_start & ~has_end).sum()
            both_mismatch = neither  # BOTH_MISMATCH is when neither matches

            samples[sample_name]['tool_disagreements'][tool]['EXACT_MATCH'] += exact
            samples[sample_name]['tool_disagreements'][tool]['START_ONLY'] += start_only
            samples[sample_name]['tool_disagreements'][tool]['END_ONLY'] += end_only
            samples[sample_name]['tool_disagreements'][tool]['BOTH_MISMATCH'] += both_mismatch

            # Pairwise
            samples[sample_name]['pairwise'][query_tool][tool]['EXACT_MATCH'] += exact
            samples[sample_name]['pairwise'][query_tool][tool]['START_ONLY'] += start_only
            samples[sample_name]['pairwise'][query_tool][tool]['END_ONLY'] += end_only
            samples[sample_name]['pairwise'][query_tool][tool]['BOTH_MISMATCH'] += both_mismatch

            # Context stratification
            for ctx in df['cds_context'].dropna().unique():
                ctx_mask = df['cds_context'] == ctx
                ctx_exact = (has_start & has_end & ctx_mask).sum()
                ctx_start = (has_start & ~has_end & ctx_mask).sum()
                ctx_end = (~has_start & has_end & ctx_mask).sum()
                ctx_neither = (~has_start & ~has_end & ctx_mask).sum()

                samples[sample_name]['context_disagreements'][tool][ctx]['EXACT_MATCH'] += ctx_exact
                samples[sample_name]['context_disagreements'][tool][ctx]['START_ONLY'] += ctx_start
                samples[sample_name]['context_disagreements'][tool][ctx]['END_ONLY'] += ctx_end
                samples[sample_name]['context_disagreements'][tool][ctx]['BOTH_MISMATCH'] += ctx_neither

        # OPTIMIZATION: Sample early by consensus score, then convert to dict
        # Group by consensus score and sample
        for score in df['consensus_score'].unique():
            score_mask = df['consensus_score'] == score
            score_df = df[score_mask]

            total_count = len(score_df)

            # Sample if needed
            if len(score_df) > max_examples_per_category:
                score_df = score_df.sample(n=max_examples_per_category, random_state=42)

            # Store consensus stats
            if 'consensus_counts' not in samples[sample_name]:
                samples[sample_name]['consensus_counts'] = {}
            samples[sample_name]['consensus_counts'][int(score)] = total_count

            # Convert sampled rows to feature dicts (still using iterrows but on much smaller subset)
            for _, row in score_df.iterrows():
                feature = {
                    'chr': row['chr'],
                    'start': int(row['start_pos']),
                    'end': int(row['end_pos']),
                    'strand': str(row.get('strand', '')),
                    'query_tool': query_tool,
                    'feature_name': str(row.get('feature_name', 'N/A')),
                    'rna_biotype': str(row.get('rna_biotype', 'intergenic')),
                    'cds_context': str(row.get('cds_context', 'intergenic')),
                    'gene_name': str(row.get('gene_name', 'intergenic')),
                    'ucsc_url': str(row.get('ucsc_url', '')),
                    'consensus_score': int(score),
                    'total_tools': len(other_tools) + 1,
                    'tool_matches': {}
                }

                # Build tool_matches for sampled features only
                for tool in other_tools:
                    s = int(row.get(f'{tool}_start_matches', 0))
                    e = int(row.get(f'{tool}_end_matches', 0))
                    cls = classify_disagreement(s, e)
                    feature['tool_matches'][tool] = {
                        'start': s,
                        'end': e,
                        'class': cls
                    }

                samples[sample_name]['features_by_consensus'][score].append(feature)

    # Finalize tool lists and consolidate sampled features
    for s in samples.values():
        s['tools'] = sorted(s['tools'])

        # Features are already sampled in the main loop, just consolidate
        sampled_features = []
        consensus_stats = {}

        for consensus_score, features in s['features_by_consensus'].items():
            sampled_features.extend(features)

            # Get actual count from consensus_counts
            total_count = s.get('consensus_counts', {}).get(consensus_score, len(features))
            # Convert numpy.int64 to Python int for JSON serialization
            consensus_stats[int(consensus_score)] = {
                'total': int(total_count),
                'sampled': len(features)
            }

        s['features'] = sampled_features
        s['consensus_stats'] = consensus_stats
        s['tool_count'] = len(s['tools'])

        # Clean up temporary structures
        del s['features_by_consensus']
        if 'consensus_counts' in s:
            del s['consensus_counts']

        total_features = sum(stats['total'] for stats in consensus_stats.values())
        print(f"Sample {s['name']}: {total_features:,} total features, "
              f"{len(sampled_features):,} sampled for report", file=sys.stderr)

    return {
        'samples': sorted(samples.values(), key=lambda x: x['name'])
    }


# ---------------------------------------------------------------------
# Per-tool statistics
# ---------------------------------------------------------------------

def compute_tool_statistics(sample):
    """
    Compute basic statistics for each tool: feature count, size distribution.
    """
    tool_stats = {}

    # Group features by query tool
    features_by_tool = defaultdict(list)
    for feature in sample['features']:
        features_by_tool[feature['query_tool']].append(feature)

    for tool in sample['tools']:
        if tool not in features_by_tool:
            # This tool might be in "other_tools" but not query_tool
            # Try to count from tool_disagreements
            total_comparisons = sum(sample['tool_disagreements'].get(tool, {}).values())
            tool_stats[tool] = {
                'feature_count': 'N/A (not query tool)',
                'mean_length': None,
                'median_length': None,
                'min_length': None,
                'max_length': None,
                'total_comparisons': total_comparisons
            }
            continue

        features = features_by_tool[tool]
        lengths = [f['end'] - f['start'] for f in features]

        tool_stats[tool] = {
            'feature_count': len(features),
            'mean_length': int(sum(lengths) / len(lengths)) if lengths else 0,
            'median_length': int(sorted(lengths)[len(lengths)//2]) if lengths else 0,
            'min_length': min(lengths) if lengths else 0,
            'max_length': max(lengths) if lengths else 0,
            'total_comparisons': sum(sample['tool_disagreements'].get(tool, {}).values())
        }

    return tool_stats


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
# Locus-level analysis
# ---------------------------------------------------------------------

def analyze_loci(sample, max_loci=50):
    """
    Group features by genomic loci and analyze agreement patterns at each locus.
    A locus is defined as overlapping features from different tools.

    Returns top N loci with interesting agreement/disagreement patterns.
    """
    from collections import namedtuple

    Interval = namedtuple('Interval', ['chr', 'start', 'end', 'strand', 'tool', 'feature'])

    # Collect all features from all tools with their coordinates
    intervals = []
    for feature in sample['features']:
        intervals.append(Interval(
            chr=feature['chr'],
            start=feature['start'],
            end=feature['end'],
            strand=feature['strand'],
            tool=feature['query_tool'],
            feature=feature
        ))

    # Sort by chromosome and position
    intervals.sort(key=lambda x: (x.chr, x.start, x.end))

    # Group overlapping intervals into loci
    loci = []
    current_locus = []

    for interval in intervals:
        if not current_locus:
            current_locus = [interval]
            continue

        # Check if this interval overlaps with current locus
        locus_chr = current_locus[0].chr
        locus_strand = current_locus[0].strand
        locus_start = min(iv.start for iv in current_locus)
        locus_end = max(iv.end for iv in current_locus)

        # Overlaps if same chr/strand and coordinates overlap
        if (interval.chr == locus_chr and
            interval.strand == locus_strand and
            interval.start < locus_end and
            interval.end > locus_start):
            current_locus.append(interval)
        else:
            # Save current locus and start new one
            if len(current_locus) > 1:  # Only keep multi-tool loci
                loci.append(current_locus)
            current_locus = [interval]

    # Don't forget last locus
    if len(current_locus) > 1:
        loci.append(current_locus)

    # Analyze each locus
    locus_analyses = []
    for locus in loci:
        tools_present = set(iv.tool for iv in locus)

        # Calculate locus boundaries
        locus_chr = locus[0].chr
        locus_strand = locus[0].strand
        locus_start = min(iv.start for iv in locus)
        locus_end = max(iv.end for iv in locus)

        # Count agreement patterns
        starts = defaultdict(list)  # start_pos -> [tools]
        ends = defaultdict(list)    # end_pos -> [tools]

        for iv in locus:
            starts[iv.start].append(iv.tool)
            ends[iv.end].append(iv.tool)

        # Determine agreement pattern
        start_consensus = max(len(tools) for tools in starts.values())
        end_consensus = max(len(tools) for tools in ends.values())
        total_tools = len(tools_present)

        # Get gene context from first feature
        gene_name = locus[0].feature.get('gene_name', 'intergenic')
        cds_context = locus[0].feature.get('cds_context', 'intergenic')
        ucsc_url = locus[0].feature.get('ucsc_url', '')

        locus_analyses.append({
            'chr': locus_chr,
            'start': locus_start,
            'end': locus_end,
            'strand': locus_strand,
            'length': locus_end - locus_start,
            'tools_present': sorted(tools_present),
            'num_tools': total_tools,
            'num_features': len(locus),
            'start_consensus': start_consensus,
            'end_consensus': end_consensus,
            'agreement_pattern': f"Start:{start_consensus}/{total_tools}, End:{end_consensus}/{total_tools}",
            'gene_name': gene_name,
            'cds_context': cds_context,
            'ucsc_url': ucsc_url,
            'features_per_tool': {tool: sum(1 for iv in locus if iv.tool == tool)
                                  for tool in tools_present}
        })

    # Sort by interesting patterns: highest tool count, then most disagreement
    locus_analyses.sort(key=lambda x: (
        -x['num_tools'],  # More tools = more interesting
        -(x['num_features'] - x['num_tools']),  # More features per tool = more disagreement
        -abs(x['start_consensus'] - x['end_consensus'])  # Start/end disagreement
    ))

    return locus_analyses[:max_loci]


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
    parser.add_argument('-m', '--max-examples', type=int, default=100,
                        help='Maximum number of feature examples per consensus level (default: 100)')

    args = parser.parse_args()

    parsed = parse_results_directory(Path(args.results_dir), max_examples_per_category=args.max_examples)
    samples = parsed['samples']

    for s in samples:
        s['tool_statistics'] = compute_tool_statistics(s)
        s['tool_signatures'] = compute_tool_signatures(s)
        s['locus_analyses'] = analyze_loci(s, max_loci=50)
        print(f"  → Found {len(s['locus_analyses'])} interesting multi-tool loci", file=sys.stderr)

    total_features_sampled = sum(len(s['features']) for s in samples)
    total_features_actual = sum(
        sum(stats['total'] for stats in s.get('consensus_stats', {}).values())
        for s in samples
    )

    # Convert all numpy types to Python native types for JSON serialization
    samples_converted = convert_numpy_types(samples)

    context = {
        'run_name': args.run_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'ucsc_session_url': args.ucsc_session_url,
        'samples': samples_converted,
        'total_samples': len(samples),
        'total_features': total_features_sampled,
        'total_features_actual': total_features_actual,
        'max_examples': args.max_examples
    }

    template = load_template(args.template)
    html = template.render(**context)

    output_path = Path(args.output)
    output_path.write_text(html)

    print(f"\n✅ HTML report generated: {output_path}", file=sys.stderr)
    print(f"   Total features analyzed: {total_features_actual:,}", file=sys.stderr)
    print(f"   Features in report: {total_features_sampled:,} (sampled)", file=sys.stderr)


if __name__ == '__main__':
    main()
