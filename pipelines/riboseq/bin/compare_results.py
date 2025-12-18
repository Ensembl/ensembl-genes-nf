#!/usr/bin/env python3
"""
Compare Ensembl riboseq results with collaborator data.

This script compares the adapter tables and frame tables from both sources
and generates a comparison report highlighting differences.

Usage:
    python compare_results.py \
        --ensembl-adapter Ensembl_adapter_results.csv \
        --ensembl-frame Ensembl_frame_results.csv \
        --collab-adapter adapter_table_meta.csv \
        --collab-frame frame_table.csv \
        --output comparison_report.txt
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def compare_adapter_tables(ensembl_df, collab_df):
    """Compare adapter-level metadata between Ensembl and collaborator data."""

    print("\n" + "="*80)
    print("ADAPTER TABLE COMPARISON")
    print("="*80)

    # Get common GSMs
    ensembl_gsms = set(ensembl_df['GSM'].unique())
    collab_gsms = set(collab_df['GSM'].unique())
    common_gsms = ensembl_gsms & collab_gsms

    print(f"\nGSM Coverage:")
    print(f"  Ensembl: {len(ensembl_gsms)} GSMs")
    print(f"  Collaborator: {len(collab_gsms)} GSMs")
    print(f"  Common: {len(common_gsms)} GSMs")
    print(f"  Ensembl only: {len(ensembl_gsms - collab_gsms)}")
    print(f"  Collaborator only: {len(collab_gsms - ensembl_gsms)}")

    if len(ensembl_gsms - collab_gsms) > 0:
        print(f"\n  Ensembl-only GSMs: {sorted(ensembl_gsms - collab_gsms)}")
    if len(collab_gsms - ensembl_gsms) > 0:
        print(f"  Collaborator-only GSMs: {sorted(collab_gsms - ensembl_gsms)}")

    # Merge on GSM for comparison
    comparison = pd.merge(
        ensembl_df, collab_df,
        on='GSM',
        suffixes=('_ensembl', '_collab'),
        how='inner'
    )

    if len(comparison) == 0:
        print("\nNo common GSMs to compare!")
        return

    # Compare numeric fields
    numeric_fields = [
        'raw_reads', 'adapter_trim', 'pct_trimmed',
        'read_length_filter', 'pct_length_filter',
        'ncRNA_contamination_reads', 'pct_ncRNA_contamination',
        'mapped_to_genome_unique', 'pct_genome_unique'
    ]

    print(f"\n\nNUMERIC FIELD COMPARISON (n={len(comparison)} GSMs)")
    print("-" * 80)

    for field in numeric_fields:
        ensembl_col = f"{field}_ensembl"
        collab_col = f"{field}_collab"

        if ensembl_col not in comparison.columns or collab_col not in comparison.columns:
            continue

        # Convert to numeric, handling empty strings
        ensembl_vals = pd.to_numeric(comparison[ensembl_col], errors='coerce')
        collab_vals = pd.to_numeric(comparison[collab_col], errors='coerce')

        # Calculate differences
        valid_mask = ~(ensembl_vals.isna() | collab_vals.isna())
        if valid_mask.sum() == 0:
            print(f"\n{field}: No valid comparisons")
            continue

        ensembl_clean = ensembl_vals[valid_mask]
        collab_clean = collab_vals[valid_mask]

        # Calculate relative difference (%)
        rel_diff = ((ensembl_clean - collab_clean) / collab_clean * 100).abs()

        print(f"\n{field}:")
        print(f"  Mean Ensembl:       {ensembl_clean.mean():>15,.2f}")
        print(f"  Mean Collaborator:  {collab_clean.mean():>15,.2f}")
        print(f"  Mean abs diff:      {(ensembl_clean - collab_clean).abs().mean():>15,.2f}")
        print(f"  Mean rel diff (%):  {rel_diff.mean():>15,.2f}%")
        print(f"  Max rel diff (%):   {rel_diff.max():>15,.2f}%")

        # Flag samples with >10% difference
        large_diff = rel_diff > 10
        if large_diff.sum() > 0:
            print(f"  Samples with >10% diff: {large_diff.sum()} / {len(rel_diff)}")
            if large_diff.sum() <= 5:
                for idx in rel_diff[large_diff].index:
                    gsm = comparison.loc[idx, 'GSM']
                    ens_val = ensembl_clean.loc[idx]
                    col_val = collab_clean.loc[idx]
                    diff_pct = rel_diff.loc[idx]
                    print(f"    {gsm}: {ens_val:.0f} vs {col_val:.0f} ({diff_pct:.1f}% diff)")

    # Compare adapter sequences
    print(f"\n\nADAPTER SEQUENCE COMPARISON")
    print("-" * 80)

    ensembl_adapters = comparison['adapter_ensembl'].fillna('')
    collab_adapters = comparison['adapter_collab'].fillna('')

    # Check if adapters match (considering truncation)
    matches = 0
    partial_matches = 0
    mismatches = 0

    for idx in comparison.index:
        ens = ensembl_adapters.loc[idx]
        col = collab_adapters.loc[idx]

        if ens == col:
            matches += 1
        elif ens and col and (ens in col or col in ens):
            partial_matches += 1
        else:
            mismatches += 1

    print(f"  Exact matches:     {matches} / {len(comparison)}")
    print(f"  Partial matches:   {partial_matches} / {len(comparison)} (one is substring of other)")
    print(f"  Mismatches:        {mismatches} / {len(comparison)}")

    if mismatches > 0 and mismatches <= 5:
        print("\n  Mismatched adapters:")
        for idx in comparison.index:
            ens = ensembl_adapters.loc[idx]
            col = collab_adapters.loc[idx]
            if ens != col and not (ens in col or col in ens):
                gsm = comparison.loc[idx, 'GSM']
                print(f"    {gsm}: '{ens}' vs '{col}'")


def compare_frame_tables(ensembl_df, collab_df):
    """Compare P-site offset frame data between Ensembl and collaborator."""

    print("\n\n" + "="*80)
    print("FRAME TABLE COMPARISON")
    print("="*80)

    # Get common GSM/read_length combinations
    ensembl_df['key'] = ensembl_df['GSM'].astype(str) + '_' + ensembl_df['read_length'].astype(str)
    collab_df['key'] = collab_df['GSM'].astype(str) + '_' + collab_df['read_length'].astype(str)

    ensembl_keys = set(ensembl_df['key'].unique())
    collab_keys = set(collab_df['key'].unique())
    common_keys = ensembl_keys & collab_keys

    print(f"\nGSM/read_length Coverage:")
    print(f"  Ensembl: {len(ensembl_keys)} combinations")
    print(f"  Collaborator: {len(collab_keys)} combinations")
    print(f"  Common: {len(common_keys)} combinations")

    # Merge for comparison
    comparison = pd.merge(
        ensembl_df, collab_df,
        on=['GSM', 'read_length'],
        suffixes=('_ensembl', '_collab'),
        how='inner'
    )

    if len(comparison) == 0:
        print("\nNo common GSM/read_length combinations to compare!")
        return

    print(f"\n\nP-SITE OFFSET COMPARISON (n={len(comparison)} combinations)")
    print("-" * 80)

    # Compare corrected_offset_from_5
    # Check which columns exist after merge
    offset_col_ens = 'corrected_offset_from_5_ensembl' if 'corrected_offset_from_5_ensembl' in comparison.columns else 'corrected_offset_from_5'
    offset_col_col = 'corrected_offset_from_5_collab' if 'corrected_offset_from_5_collab' in comparison.columns else 'p_site_shift'

    ensembl_offsets = pd.to_numeric(comparison[offset_col_ens], errors='coerce')
    collab_offsets = pd.to_numeric(comparison[offset_col_col], errors='coerce')

    valid_mask = ~(ensembl_offsets.isna() | collab_offsets.isna())
    ensembl_clean = ensembl_offsets[valid_mask]
    collab_clean = collab_offsets[valid_mask]

    # Offset differences
    offset_diff = (ensembl_clean - collab_clean).abs()
    matches = (offset_diff == 0).sum()

    print(f"\ncorrected_offset_from_5:")
    print(f"  Exact matches:     {matches} / {len(ensembl_clean)} ({matches/len(ensembl_clean)*100:.1f}%)")
    print(f"  Mean abs diff:     {offset_diff.mean():.2f}")
    print(f"  Max abs diff:      {offset_diff.max():.0f}")

    # Show distribution of differences
    diff_counts = offset_diff.value_counts().sort_index()
    print(f"\n  Offset difference distribution:")
    for diff, count in diff_counts.head(10).items():
        print(f"    Diff {diff:.0f}: {count} cases ({count/len(offset_diff)*100:.1f}%)")

    # Compare total_percentage
    print(f"\n\ntotal_percentage:")
    pct_col_ens = 'total_percentage_ensembl' if 'total_percentage_ensembl' in comparison.columns else 'total_percentage'
    pct_col_col = 'total_percentage_collab' if 'total_percentage_collab' in comparison.columns else 'total_percentage'

    ensembl_pct = pd.to_numeric(comparison[pct_col_ens], errors='coerce')
    collab_pct = pd.to_numeric(comparison[pct_col_col], errors='coerce')

    valid_mask = ~(ensembl_pct.isna() | collab_pct.isna())
    if valid_mask.sum() > 0:
        ensembl_pct_clean = ensembl_pct[valid_mask]
        collab_pct_clean = collab_pct[valid_mask]

        rel_diff = ((ensembl_pct_clean - collab_pct_clean) / collab_pct_clean * 100).abs()

        print(f"  Mean Ensembl:      {ensembl_pct_clean.mean():.3f}%")
        print(f"  Mean Collaborator: {collab_pct_clean.mean():.3f}%")
        print(f"  Mean rel diff:     {rel_diff.mean():.2f}%")
        print(f"  Max rel diff:      {rel_diff.max():.2f}%")

    # Breakdown by read length
    print(f"\n\nOFFSET AGREEMENT BY READ LENGTH:")
    print("-" * 80)

    for read_len in sorted(comparison['read_length'].unique()):
        subset = comparison[comparison['read_length'] == read_len]

        ens_off = pd.to_numeric(subset[offset_col_ens], errors='coerce')
        col_off = pd.to_numeric(subset[offset_col_col], errors='coerce')

        valid = ~(ens_off.isna() | col_off.isna())
        if valid.sum() == 0:
            continue

        matches = ((ens_off[valid] - col_off[valid]).abs() == 0).sum()
        total = valid.sum()

        print(f"  Length {read_len:2d}: {matches:3d} / {total:3d} matches ({matches/total*100:5.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description='Compare Ensembl riboseq results with collaborator data'
    )
    parser.add_argument(
        '--ensembl-adapter',
        type=Path,
        required=True,
        help='Ensembl adapter table CSV'
    )
    parser.add_argument(
        '--ensembl-frame',
        type=Path,
        required=True,
        help='Ensembl frame table CSV'
    )
    parser.add_argument(
        '--collab-adapter',
        type=Path,
        required=True,
        help='Collaborator adapter table CSV'
    )
    parser.add_argument(
        '--collab-frame',
        type=Path,
        required=True,
        help='Collaborator frame table CSV'
    )
    parser.add_argument(
        '--output',
        type=Path,
        help='Output report file (default: print to stdout)'
    )

    args = parser.parse_args()

    # Read data
    print("Loading data...")
    ensembl_adapter = pd.read_csv(args.ensembl_adapter)
    ensembl_frame = pd.read_csv(args.ensembl_frame)
    collab_adapter = pd.read_csv(args.collab_adapter)
    collab_frame = pd.read_csv(args.collab_frame)

    print(f"  Ensembl adapter: {len(ensembl_adapter)} rows")
    print(f"  Ensembl frame:   {len(ensembl_frame)} rows")
    print(f"  Collab adapter:  {len(collab_adapter)} rows")
    print(f"  Collab frame:    {len(collab_frame)} rows")

    # Compare
    compare_adapter_tables(ensembl_adapter, collab_adapter)
    compare_frame_tables(ensembl_frame, collab_frame)

    print("\n" + "="*80)
    print("COMPARISON COMPLETE")
    print("="*80 + "\n")


if __name__ == '__main__':
    main()
