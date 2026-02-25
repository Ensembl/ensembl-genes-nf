#!/usr/bin/env python3
"""
Collect riboseq workflow output data and format it to match collaborator's CSV format.

This script parses the published outputs from the riboseq workflow and generates
two CSV files:

1. Adapter Table (sample-level metadata):
   - Specie
   - GSM
   - GSE
   - sample_name
   - raw_reads
   - adapter (from adapter_report.fa or fastp JSON)
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

2. Frame Table (per-read-length P-site offset data):
   - GSE
   - GSM
   - sample_name
   - mapped_to_genome_unique
   - read_length
   - total_percentage
   - corrected_offset_from_5

Usage:
    # Basic usage (no aggregation)
    python collect_comparison_data.py --outdir /path/to/workflow/results

    # Custom output filenames
    python collect_comparison_data.py --outdir /path/to/workflow/results \
        --adapter-table my_adapter.csv --frame-table my_frame.csv

    # Aggregate multiple SRR runs to GSM level
    python collect_comparison_data.py --outdir /path/to/workflow/results \
        --gsm-mapping gsm_mapping.csv

Aggregation Strategy (when --gsm-mapping is provided):
    - Count fields (raw_reads, adapter_trim, etc.): sum across runs
    - Percentage fields: recalculated from summed counts
    - Adapter sequences: most common adapter
    - P-site offsets: offset from run with highest read abundance per length
"""

import argparse
import json
import csv
import re
import gzip
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


def parse_adapter_report(fa_path: Path) -> str:
    """Parse adapter report FASTA file to extract detected adapter sequence.

    Args:
        fa_path: Path to the adapter_report.fa file

    Returns:
        Adapter sequence string
    """
    adapter_seq = ""

    try:
        with open(fa_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    # This is the sequence line
                    adapter_seq = line
                    break
    except Exception as e:
        print(f"Warning: Error parsing {fa_path}: {e}")

    return adapter_seq


def parse_ribowaltz_psite_offset(tsv_path: Path) -> Dict:
    """Parse ribowaltz P-site offset TSV file to extract offset information.

    Args:
        tsv_path: Path to the .psite_offset.tsv or .psite_offset.tsv.gz file

    Returns:
        Dictionary mapping read length to offset information
    """
    result = {}

    # Handle both gzipped and uncompressed files
    if tsv_path.suffix == '.gz':
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    try:
        with open_func(tsv_path, mode) as f:
            # Skip header
            _ = f.readline()

            # Parse each line
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                length = fields[0]
                total_pct = fields[1]
                start_pct = fields[2]
                around_start = fields[3]
                offset_from_5 = fields[4]
                offset_from_3 = fields[5]
                corrected_offset_from_5 = fields[6]
                corrected_offset_from_3 = fields[7]

                result[length] = {
                    'total_percentage': total_pct,
                    'start_percentage': start_pct,
                    'around_start': around_start,
                    'offset_from_5': offset_from_5,
                    'offset_from_3': offset_from_3,
                    'corrected_offset_from_5': corrected_offset_from_5,
                    'corrected_offset_from_3': corrected_offset_from_3
                }
    except Exception as e:
        print(f"Warning: Error parsing {tsv_path}: {e}")

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
            'UMI_seq_3': '',
            'psite_offsets': ''  # Will store offsets as JSON or formatted string
        }

        # Parse adapter report (if available) to get detected adapter
        adapter_report_dir = outdir / "adapter_reports"
        adapter_report = adapter_report_dir / f"{sample_id}_adapter_report.fa"
        if not adapter_report.exists():
            # Try alternative naming
            adapter_reports = list(adapter_report_dir.glob(f"*{sample_id}*_adapter_report.fa"))
            if adapter_reports:
                adapter_report = adapter_reports[0]

        if adapter_report.exists():
            try:
                detected_adapter = parse_adapter_report(adapter_report)
                if detected_adapter:
                    row['adapter'] = detected_adapter
            except Exception as e:
                print(f"Warning: Error parsing {adapter_report}: {e}")

        # Parse provided adapter fastp JSON
        try:
            fastp_provided = parse_fastp_json(provided_json)
            row['raw_reads'] = fastp_provided['total_reads_before']
            # Use fastp adapter if we don't have one from adapter report
            if not row['adapter']:
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

        # Parse ribowaltz P-site offset data
        ribowaltz_dir = outdir / "ribowaltz"
        psite_offset_file = ribowaltz_dir / f"{sample_id}.psite_offset.tsv.gz"
        if not psite_offset_file.exists():
            # Try uncompressed version
            psite_offset_file = ribowaltz_dir / f"{sample_id}.psite_offset.tsv"
        if not psite_offset_file.exists():
            # Try alternative naming patterns
            psite_files = list(ribowaltz_dir.glob(f"*{sample_id}*.psite_offset.tsv*"))
            if psite_files:
                psite_offset_file = psite_files[0]

        if psite_offset_file.exists():
            try:
                psite_data = parse_ribowaltz_psite_offset(psite_offset_file)
                # Format as JSON string for storage
                row['psite_offsets'] = json.dumps(psite_data)
            except Exception as e:
                print(f"Warning: Error parsing {psite_offset_file}: {e}")

        results.append(row)

    return results


def write_adapter_table(data: List[Dict], output_path: Path):
    """Write sample-level metadata to CSV file (adapter table format)."""
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

    # Write data without psite_offsets (sample-level only)
    rows_to_write = []
    for row in data:
        row_copy = {k: v for k, v in row.items() if k in fieldnames}
        rows_to_write.append(row_copy)

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_to_write)

    print(f"Wrote {len(rows_to_write)} rows to {output_path}")


def aggregate_runs_to_gsm(data: List[Dict], gsm_mapping: Dict[str, str]) -> List[Dict]:
    """
    Aggregate multiple SRR runs to their parent GSM accessions.

    For count fields: sum across runs
    For percentage fields: recalculate from summed counts
    For adapter sequences: use most common
    For P-site offsets: use most abundant read length's offset

    Args:
        data: List of sample dictionaries (one per SRR)
        gsm_mapping: Dict mapping sample_id (SRR) to GSM

    Returns:
        List of aggregated dictionaries (one per GSM)
    """
    from collections import defaultdict

    # Group runs by GSM
    gsm_groups = defaultdict(list)
    for row in data:
        sample_id = row.get('sample_name', '')
        gsm = gsm_mapping.get(sample_id, row.get('GSM', ''))
        gsm_groups[gsm].append(row)

    aggregated = []

    for gsm, runs in gsm_groups.items():
        if len(runs) == 1:
            # Single run, no aggregation needed
            aggregated.append(runs[0])
            continue

        # Multiple runs - aggregate
        agg_row = {
            'Specie': runs[0].get('Specie', ''),
            'GSM': gsm,
            'GSE': runs[0].get('GSE', ''),
            'sample_name': f"{gsm} (n={len(runs)} runs)",
            'UMI': runs[0].get('UMI', 'No'),
            'UMI_seq_5': runs[0].get('UMI_seq_5', ''),
            'UMI_seq_3': runs[0].get('UMI_seq_3', '')
        }

        # Count fields - sum across runs (absolute counts)
        count_fields = [
            'raw_reads', 'adapter_trim', 'read_length_filter',
            'ncRNA_contamination_reads', 'mapped_to_genome_unique'
        ]

        for field in count_fields:
            values = [r.get(field, 0) for r in runs if r.get(field) not in ['', None]]
            if values:
                # Convert to float and sum
                numeric_values = [float(v) if v != '' else 0 for v in values]
                agg_row[field] = sum(numeric_values)
            else:
                agg_row[field] = ''

        # Percentage fields - recalculate from summed counts
        if agg_row['raw_reads'] and agg_row['adapter_trim']:
            agg_row['pct_trimmed'] = agg_row['adapter_trim'] / agg_row['raw_reads']
        else:
            agg_row['pct_trimmed'] = ''

        if agg_row.get('adapter_trim') and agg_row.get('read_length_filter'):
            agg_row['pct_length_filter'] = agg_row['read_length_filter'] / agg_row['adapter_trim']
        else:
            agg_row['pct_length_filter'] = ''

        if agg_row.get('read_length_filter') and agg_row.get('ncRNA_contamination_reads'):
            agg_row['pct_ncRNA_contamination'] = (agg_row['ncRNA_contamination_reads'] / agg_row['read_length_filter']) * 100
        else:
            agg_row['pct_ncRNA_contamination'] = ''

        # Calculate pct_genome_unique from summed values
        if agg_row.get('read_length_filter') and agg_row.get('ncRNA_contamination_reads') and agg_row.get('mapped_to_genome_unique'):
            post_rrna_reads = agg_row['read_length_filter'] - agg_row['ncRNA_contamination_reads']
            if post_rrna_reads > 0:
                agg_row['pct_genome_unique'] = (agg_row['mapped_to_genome_unique'] / post_rrna_reads) * 100
            else:
                agg_row['pct_genome_unique'] = ''
        elif agg_row.get('read_length_filter') and agg_row.get('mapped_to_genome_unique'):
            agg_row['pct_genome_unique'] = (agg_row['mapped_to_genome_unique'] / agg_row['read_length_filter']) * 100
        else:
            agg_row['pct_genome_unique'] = ''

        # Adapter - use most common
        adapters = [r.get('adapter', '') for r in runs if r.get('adapter')]
        if adapters:
            from collections import Counter
            agg_row['adapter'] = Counter(adapters).most_common(1)[0][0]
        else:
            agg_row['adapter'] = ''

        # P-site offsets - merge all offsets, weighted by total_percentage
        merged_offsets = {}
        for run in runs:
            if not run.get('psite_offsets'):
                continue
            try:
                offsets = json.loads(run['psite_offsets'])
                for length, offset_data in offsets.items():
                    if length not in merged_offsets:
                        merged_offsets[length] = []
                    merged_offsets[length].append(offset_data)
            except (json.JSONDecodeError, TypeError):
                continue

        # For each read length, use the offset from the run with highest total_percentage
        final_offsets = {}
        for length, offset_list in merged_offsets.items():
            # Find offset with max total_percentage
            best_offset = max(offset_list,
                            key=lambda x: float(x.get('total_percentage', 0)) if x.get('total_percentage') else 0)
            final_offsets[length] = best_offset

        agg_row['psite_offsets'] = json.dumps(final_offsets) if final_offsets else ''

        aggregated.append(agg_row)

    return aggregated


def write_frame_table(data: List[Dict], output_path: Path):
    """Write per-read-length P-site offset data to CSV file (frame table format)."""
    if not data:
        print("Warning: No data to write")
        return

    fieldnames = [
        'GSE', 'GSM', 'sample_name', 'mapped_to_genome_unique',
        'read_length', 'total_percentage', 'corrected_offset_from_5'
    ]

    rows_to_write = []
    for row in data:
        if not row.get('psite_offsets'):
            continue

        # Parse the JSON offset data
        try:
            offsets = json.loads(row['psite_offsets'])
        except (json.JSONDecodeError, TypeError):
            continue

        # Create one row per read length
        for length, offset_data in offsets.items():
            frame_row = {
                'GSE': row.get('GSE', ''),
                'GSM': row.get('GSM', ''),
                'sample_name': row.get('sample_name', ''),
                'mapped_to_genome_unique': row.get('mapped_to_genome_unique', ''),
                'read_length': length,
                'total_percentage': offset_data.get('total_percentage', ''),
                'corrected_offset_from_5': offset_data.get('corrected_offset_from_5', '')
            }
            rows_to_write.append(frame_row)

    if rows_to_write:
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows_to_write)
        print(f"Wrote {len(rows_to_write)} rows to {output_path}")
    else:
        print(f"Warning: No P-site offset data found, skipping {output_path}")


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
        '--adapter-table',
        type=Path,
        default='adapter_table.csv',
        help='Output CSV file for sample-level metadata (default: adapter_table.csv)'
    )
    parser.add_argument(
        '--frame-table',
        type=Path,
        default='frame_table.csv',
        help='Output CSV file for per-read-length P-site offsets (default: frame_table.csv)'
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
    parser.add_argument(
        '--gsm-mapping',
        type=Path,
        help='Optional CSV file mapping SRR to GSM (columns: Run, GSM). If provided, multiple runs will be aggregated to GSM level.'
    )

    args = parser.parse_args()

    # Validate output directory
    if not args.outdir.exists():
        parser.error(f"Output directory not found: {args.outdir}")

    # Collect data
    print(f"Collecting data from {args.outdir}...")
    data = collect_run_data(args.outdir, species=args.species, gse=args.gse)

    # If GSM mapping provided, aggregate runs to GSM level
    if args.gsm_mapping:
        if not args.gsm_mapping.exists():
            parser.error(f"GSM mapping file not found: {args.gsm_mapping}")

        print(f"Loading GSM mapping from {args.gsm_mapping}...")
        import pandas as pd
        mapping_df = pd.read_csv(args.gsm_mapping)

        # Create dict mapping Run -> GSM
        gsm_mapping = dict(zip(mapping_df['Run'], mapping_df['GSM']))
        print(f"Found {len(gsm_mapping)} Run->GSM mappings")

        print(f"Aggregating {len(data)} runs to GSM level...")
        data = aggregate_runs_to_gsm(data, gsm_mapping)
        print(f"After aggregation: {len(data)} samples")

    # Write both output files
    write_adapter_table(data, args.adapter_table)
    write_frame_table(data, args.frame_table)


if __name__ == '__main__':
    main()
