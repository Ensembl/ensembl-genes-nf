#!/usr/bin/env python3
"""
Collect riboseq workflow output data and format it to match collaborator's CSV format.

This script parses the published outputs from the riboseq workflow and generates
a CSV file with the following columns:
- Specie
- GSM
- GSE
- sample_name
- raw_reads
- adapter
- adapter_trim
- pct_trimmed
- read_length_filter
- pct_length_filter
- ncRNA_contamination_reads
- pct_ncRNA_contamination
- mapped_to_genome_unique
- pct_genome_unique
- UMI
- UMI_seq_5
- UMI_seq_3

Usage:
    python collect_comparison_data.py --outdir /path/to/workflow/results --output comparison.csv
"""

import argparse
import json
import csv
import re
from pathlib import Path
from typing import Dict, Optional, List


def parse_fastp_json(json_path: Path) -> Dict:
    """Parse fastp JSON output to extract trimming statistics."""
    with open(json_path, 'r') as f:
        data = json.load(f)

    result = {}

    # Summary stats
    summary = data.get('summary', {})
    before = summary.get('before_filtering', {})
    after = summary.get('after_filtering', {})

    result['total_reads_before'] = before.get('total_reads', 0)
    result['total_reads_after'] = after.get('total_reads', 0)

    # Adapter trimming
    adapter_cutting = data.get('adapter_cutting', {})
    result['adapter_trimmed_reads'] = adapter_cutting.get('adapter_trimmed_reads', 0)
    result['adapter_sequence'] = adapter_cutting.get('adapter_sequence', '')

    # Calculate percentages
    if result['total_reads_before'] > 0:
        result['pct_trimmed'] = result['adapter_trimmed_reads'] / result['total_reads_before']
        result['pct_length_filter'] = result['total_reads_after'] / result['adapter_trimmed_reads'] if result['adapter_trimmed_reads'] > 0 else 0
    else:
        result['pct_trimmed'] = 0
        result['pct_length_filter'] = 0

    return result


def parse_rrna_filter_log(log_path: Path) -> Dict:
    """Parse rRNA filter log to extract contamination statistics."""
    result = {}

    with open(log_path, 'r') as f:
        content = f.read()

    # Extract rRNA reads
    match = re.search(r'rRNA reads:\s+(\d+)', content)
    if match:
        result['rrna_reads'] = int(match.group(1))

    # Extract percentage rRNA
    match = re.search(r'Percentage rRNA:\s+([\d.]+)', content)
    if match:
        result['pct_rrna'] = float(match.group(1))

    # Extract filtered reads
    match = re.search(r'Filtered reads \(no rRNA\):\s+(\d+)', content)
    if match:
        result['filtered_reads'] = int(match.group(1))

    return result


def parse_star_log(log_path: Path) -> Dict:
    """Parse STAR Log.final.out to extract alignment statistics."""
    result = {}

    with open(log_path, 'r') as f:
        content = f.read()

    # Extract number of input reads
    match = re.search(r'Number of input reads \|\s+(\d+)', content)
    if match:
        result['input_reads'] = int(match.group(1))

    # Extract uniquely mapped reads
    match = re.search(r'Uniquely mapped reads number \|\s+(\d+)', content)
    if match:
        result['uniquely_mapped'] = int(match.group(1))

    # Extract uniquely mapped percentage
    match = re.search(r'Uniquely mapped reads % \|\s+([\d.]+)%', content)
    if match:
        result['uniquely_mapped_pct'] = float(match.group(1))

    # Extract multi-mapped reads
    match = re.search(r'Number of reads mapped to multiple loci \|\s+(\d+)', content)
    if match:
        result['multi_mapped'] = int(match.group(1))

    # Extract multi-mapped percentage
    match = re.search(r'% of reads mapped to multiple loci \|\s+([\d.]+)%', content)
    if match:
        result['multi_mapped_pct'] = float(match.group(1))

    return result


def extract_sample_id(filename: str) -> str:
    """Extract sample ID from filename (typically SRR* accession)."""
    # Common patterns: SRR14890797, GSM1299091, etc.
    match = re.search(r'(SRR\d+|ERR\d+|DRR\d+|GSM\d+)', filename)
    if match:
        return match.group(1)
    # Fallback: use filename stem
    return Path(filename).stem.split('_')[0]


def collect_run_data(outdir: Path, species: str = "", gse: str = "") -> List[Dict]:
    """
    Collect data for all runs in the output directory.

    Args:
        outdir: Path to workflow output directory
        species: Species name (default: blank)
        gse: GSE accession (default: blank)

    Returns:
        List of dictionaries, one per run
    """
    results = []

    # Find all fastp JSON files (provided adapter version)
    fastp_dir = outdir / "fastp"
    if not fastp_dir.exists():
        print(f"Warning: {fastp_dir} not found")
        return results

    # Process each sample
    for provided_json in sorted(fastp_dir.glob("*_provided_fastp.json")):
        sample_id = extract_sample_id(provided_json.name)

        row = {
            'Specie': species,
            'GSM': sample_id,  # Or blank if not available
            'GSE': gse,
            'sample_name': sample_id,  # Could be customized
            'raw_reads': '',
            'adapter': '',
            'adapter_trim': '',
            'pct_trimmed': '',
            'read_length_filter': '',
            'pct_length_filter': '',
            'ncRNA_contamination_reads': '',
            'pct_ncRNA_contamination': '',
            'mapped_to_genome_unique': '',
            'pct_genome_unique': '',
            'UMI': 'No',
            'UMI_seq_5': '',
            'UMI_seq_3': ''
        }

        # Parse provided adapter fastp JSON
        try:
            fastp_provided = parse_fastp_json(provided_json)
            row['raw_reads'] = fastp_provided['total_reads_before']
            row['adapter'] = fastp_provided['adapter_sequence']
            row['adapter_trim'] = fastp_provided['adapter_trimmed_reads']
            row['pct_trimmed'] = fastp_provided['pct_trimmed']

            # Parse auto-detected adapter fastp JSON (for length filtering)
            auto_json = provided_json.parent / provided_json.name.replace('_provided_fastp.json', '_auto_fastp.json')
            if auto_json.exists():
                fastp_auto = parse_fastp_json(auto_json)
                row['read_length_filter'] = fastp_auto['total_reads_after']
                # Calculate percentage of reads passing length filter relative to trimmed reads
                if row['adapter_trim'] > 0:
                    row['pct_length_filter'] = fastp_auto['total_reads_after'] / row['adapter_trim']
        except Exception as e:
            print(f"Warning: Error parsing {provided_json}: {e}")

        # Parse rRNA filter log (if available)
        rrna_log_dir = outdir / "rrna_filter"
        rrna_log = rrna_log_dir / f"{sample_id}_rrna_filter.log"
        if not rrna_log.exists():
            # Try alternative naming
            rrna_logs = list(rrna_log_dir.glob(f"*{sample_id}*_rrna_filter.log"))
            if rrna_logs:
                rrna_log = rrna_logs[0]

        if rrna_log.exists():
            try:
                rrna_data = parse_rrna_filter_log(rrna_log)
                row['ncRNA_contamination_reads'] = rrna_data.get('rrna_reads', '')
                row['pct_ncRNA_contamination'] = rrna_data.get('pct_rrna', '')
            except Exception as e:
                print(f"Warning: Error parsing {rrna_log}: {e}")

        # Parse STAR alignment log
        star_log_dir = outdir / "star_align" / "logs"
        star_log = star_log_dir / f"{sample_id}.Log.final.out"
        if not star_log.exists():
            # Try alternative naming
            star_logs = list(star_log_dir.glob(f"*{sample_id}*.Log.final.out"))
            if star_logs:
                star_log = star_logs[0]

        if star_log.exists():
            try:
                star_data = parse_star_log(star_log)
                row['mapped_to_genome_unique'] = star_data.get('uniquely_mapped', '')
                # Calculate percentage of unique mapping relative to rRNA-filtered reads (if available) or length-filtered reads
                input_for_mapping = row['ncRNA_contamination_reads']
                if input_for_mapping and row['read_length_filter']:
                    # If rRNA filtering was done, calculate based on post-rRNA reads
                    post_rrna_reads = row['read_length_filter'] - input_for_mapping
                    if post_rrna_reads > 0 and row['mapped_to_genome_unique']:
                        row['pct_genome_unique'] = (row['mapped_to_genome_unique'] / post_rrna_reads) * 100
                elif row['read_length_filter'] and row['mapped_to_genome_unique']:
                    # No rRNA filtering, calculate based on length-filtered reads
                    row['pct_genome_unique'] = (row['mapped_to_genome_unique'] / row['read_length_filter']) * 100
            except Exception as e:
                print(f"Warning: Error parsing {star_log}: {e}")

        results.append(row)

    return results


def write_csv(data: List[Dict], output_path: Path):
    """Write collected data to CSV file."""
    if not data:
        print("Warning: No data to write")
        return

    fieldnames = [
        'Specie', 'GSM', 'GSE', 'sample_name', 'raw_reads', 'adapter',
        'adapter_trim', 'pct_trimmed', 'read_length_filter', 'pct_length_filter',
        'ncRNA_contamination_reads', 'pct_ncRNA_contamination',
        'mapped_to_genome_unique', 'pct_genome_unique',
        'UMI', 'UMI_seq_5', 'UMI_seq_3'
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

    print(f"Wrote {len(data)} rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Collect riboseq workflow data and format for comparison with collaborators'
    )
    parser.add_argument(
        '--outdir',
        type=Path,
        required=True,
        help='Path to workflow output directory (contains fastp/, star_align/, etc.)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default='comparison_table.csv',
        help='Output CSV file path (default: comparison_table.csv)'
    )
    parser.add_argument(
        '--species',
        type=str,
        default='',
        help='Species name (e.g., "Zebrafish")'
    )
    parser.add_argument(
        '--gse',
        type=str,
        default='',
        help='GSE accession (e.g., "GSE53693")'
    )

    args = parser.parse_args()

    # Validate output directory
    if not args.outdir.exists():
        parser.error(f"Output directory not found: {args.outdir}")

    # Collect data
    print(f"Collecting data from {args.outdir}...")
    data = collect_run_data(args.outdir, species=args.species, gse=args.gse)

    # Write CSV
    write_csv(data, args.output)


if __name__ == '__main__':
    main()
