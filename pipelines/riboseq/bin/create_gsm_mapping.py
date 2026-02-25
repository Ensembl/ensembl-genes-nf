#!/usr/bin/env python3
"""
Create SRR to GSM mapping file using pysradb.

This script takes a list of GSM accessions and creates a mapping file
showing which SRR runs correspond to each GSM sample.

Usage:
    python create_gsm_mapping.py --input gsm_list.txt --output gsm_mapping.csv
    python create_gsm_mapping.py --gsm GSM1299091 GSM1299092 --output gsm_mapping.csv
"""

import argparse
import sys
from pathlib import Path

try:
    from pysradb import SRAweb
    import pandas as pd
except ImportError:
    print("Error: pysradb is required. Install with: pip install pysradb")
    sys.exit(1)


def create_gsm_mapping(gsm_list, output_csv='gsm_mapping.csv'):
    """
    Create SRR to GSM mapping file.

    Parameters:
    -----------
    gsm_list : list
        List of GSM accessions
    output_csv : str
        Path to output CSV file
    """
    print(f"Processing {len(gsm_list)} GSM accessions...")

    db = SRAweb()
    gsm_to_srr = db.gsm_to_srr(gsm_list)

    # Debug: print available columns
    print(f"Available columns: {list(gsm_to_srr.columns)}")

    # Map common column name variations
    possible_gsm_cols = ['experiment_accession', 'experiment_alias', 'sample_accession', 'sample_alias']
    possible_run_cols = ['run_accession', 'run']

    gsm_col = None
    run_col = None

    # Find GSM column
    for col in possible_gsm_cols:
        if col in gsm_to_srr.columns:
            gsm_col = col
            break

    # Find Run column
    for col in possible_run_cols:
        if col in gsm_to_srr.columns:
            run_col = col
            break

    # If standard names not found, try fuzzy matching
    if gsm_col is None:
        for col in gsm_to_srr.columns:
            col_lower = col.lower()
            if any(x in col_lower for x in ['experiment', 'sample', 'gsm']):
                gsm_col = col
                break

    if run_col is None:
        for col in gsm_to_srr.columns:
            col_lower = col.lower()
            if 'run' in col_lower:
                run_col = col
                break

    if gsm_col is None or run_col is None:
        print("\nError: Could not find expected columns")
        print(f"Available columns: {list(gsm_to_srr.columns)}")
        print("\nFull dataframe:")
        print(gsm_to_srr.head())
        sys.exit(1)

    print(f"Using columns: GSM={gsm_col}, Run={run_col}")

    # Create mapping with Run and GSM columns
    mapping = gsm_to_srr[[gsm_col, run_col]].copy()
    mapping.columns = ['GSM', 'Run']
    mapping = mapping.drop_duplicates().sort_values('Run')

    mapping.to_csv(output_csv, index=False)

    print(f"\nMapping saved to {output_csv}")
    print(f"Total runs: {len(mapping)}")
    print(f"Unique GSMs: {len(mapping['GSM'].unique())}")

    # Show summary
    runs_per_gsm = mapping.groupby('GSM').size()
    print(f"\nRuns per GSM:")
    print(f"  Min: {runs_per_gsm.min()}")
    print(f"  Max: {runs_per_gsm.max()}")
    print(f"  Mean: {runs_per_gsm.mean():.2f}")

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description='Create SRR to GSM mapping file using pysradb'
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--input',
        type=Path,
        help='Input file with one GSM per line'
    )
    input_group.add_argument(
        '--gsm',
        nargs='+',
        help='GSM accessions (space-separated)'
    )

    parser.add_argument(
        '--output',
        type=Path,
        default='gsm_mapping.csv',
        help='Output CSV file (default: gsm_mapping.csv)'
    )

    args = parser.parse_args()

    # Get GSM list
    if args.input:
        if not args.input.exists():
            parser.error(f"Input file not found: {args.input}")
        with open(args.input) as f:
            gsm_list = [line.strip() for line in f if line.strip()]
    else:
        gsm_list = args.gsm

    # Create mapping
    mapping = create_gsm_mapping(gsm_list, args.output)

    print("\nFirst few rows:")
    print(mapping.head(10))


if __name__ == '__main__':
    main()
