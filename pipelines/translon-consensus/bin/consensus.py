#!/usr/bin/env python3
"""
Consensus analysis for ncORF annotations across multiple tools.
Computes start/end/both match counts for each feature against all other tools.
"""

import argparse
import sys
from pathlib import Path
import duckdb
import pandas as pd


def parse_samplesheet(samplesheet_path):
    """
    Parse samplesheet TSV with columns: tool, bed_path
    Returns: dict of {tool: bed_path}
    """
    df = pd.read_csv(samplesheet_path, sep='\t', header=0)
    
    if not {'tool', 'bed_path'}.issubset(df.columns):
        raise ValueError("Samplesheet must have columns: tool, bed_path")
    
    # Filter out missing files
    valid_entries = {}
    for _, row in df.iterrows():
        tool = row['tool']
        bed_path = Path(row['bed_path'])
        
        if not bed_path.exists():
            print(f"Warning: BED file not found for {tool}: {bed_path}", file=sys.stderr)
            continue
        
        valid_entries[tool] = str(bed_path)
    
    if len(valid_entries) < 2:
        raise ValueError(f"Need at least 2 valid tool BED files, found {len(valid_entries)}")
    
    return valid_entries


def ingest_beds(con, sample_name, bed_files_dict):
    """
    Load all BED files into DuckDB features table.
    
    Args:
        con: DuckDB connection
        sample_name: Sample identifier
        bed_files_dict: {tool: bed_path}
    """
    # Create table if not exists
    con.execute('''
        CREATE TABLE IF NOT EXISTS features (
            sample TEXT,
            tool TEXT,
            chr TEXT,
            start_pos INTEGER,
            end_pos INTEGER,
            strand TEXT,
            feature_name TEXT,
            score REAL,
            thickStart INTEGER,
            thickEnd INTEGER,
            itemRgb TEXT,
            blockCount INTEGER,
            blockSizes TEXT,
            blockStarts TEXT
        )
    ''')
    
    for tool, bed_path in bed_files_dict.items():
        print(f"Ingesting {tool}: {bed_path}", file=sys.stderr)
        
        # Try to infer BED format (BED3, BED6, BED12)
        # Read first line to count fields
        with open(bed_path) as f:
            first_line = f.readline().strip()
            if not first_line or first_line.startswith('#') or first_line.startswith('track'):
                # Skip header/track lines
                first_line = f.readline().strip()
            n_fields = len(first_line.split('\t'))
        
        # Build column spec based on field count
        if n_fields >= 12:
            columns = {
                'chr': 'VARCHAR', 'start_pos': 'INTEGER', 'end_pos': 'INTEGER',
                'feature_name': 'VARCHAR', 'score': 'REAL', 'strand': 'VARCHAR',
                'thickStart': 'INTEGER', 'thickEnd': 'INTEGER', 'itemRgb': 'VARCHAR',
                'blockCount': 'INTEGER', 'blockSizes': 'VARCHAR', 'blockStarts': 'VARCHAR'
            }
        elif n_fields >= 6:
            columns = {
                'chr': 'VARCHAR', 'start_pos': 'INTEGER', 'end_pos': 'INTEGER',
                'feature_name': 'VARCHAR', 'score': 'REAL', 'strand': 'VARCHAR'
            }
        elif n_fields >= 3:
            columns = {
                'chr': 'VARCHAR', 'start_pos': 'INTEGER', 'end_pos': 'INTEGER'
            }
        else:
            raise ValueError(f"Invalid BED format in {bed_path}: only {n_fields} fields")
        
        col_spec = {f'col{i}': dtype for i, dtype in enumerate(columns.values())}
        col_names = list(columns.keys())
        
        # Insert with tool and sample labels
        select_cols = ', '.join([f'col{i} as {name}' for i, name in enumerate(col_names)])
        
        try:
            con.execute(f"""
                INSERT INTO features (sample, tool, {', '.join(col_names)})
                SELECT 
                    '{sample_name}' as sample,
                    '{tool}' as tool,
                    {select_cols}
                FROM read_csv_auto('{bed_path}', 
                    delim='\t',
                    header=false,
                    columns={col_spec},
                    ignore_errors=true)
            """)
        except Exception as e:
            print(f"Error ingesting {tool} from {bed_path}: {e}", file=sys.stderr)
            raise


def compute_matches(con, sample_name, query_tool, other_tools):
    """
    For each feature in query_tool, count matches in all other tools.
    
    Args:
        con: DuckDB connection
        sample_name: Sample identifier
        query_tool: Tool to analyze
        other_tools: List of other tool names to compare against
        
    Returns:
        DataFrame with original features + match count columns
    """
    # Build lateral join clauses for each other tool
    lateral_joins = []
    select_cols = []
    
    for other_tool in other_tools:
        # Sanitize tool name for column names
        safe_name = other_tool.replace('-', '_').replace('.', '_')
        
        lateral_joins.append(f"""
            LEFT JOIN LATERAL (
                SELECT 
                    COUNT(*) FILTER (WHERE start_pos = q.start_pos) as {safe_name}_start_matches,
                    COUNT(*) FILTER (WHERE end_pos = q.end_pos) as {safe_name}_end_matches,
                    COUNT(*) FILTER (WHERE start_pos = q.start_pos AND end_pos = q.end_pos) as {safe_name}_both_matches
                FROM features
                WHERE sample = '{sample_name}' 
                  AND tool = '{other_tool}'
                  AND chr = q.chr 
                  AND strand = q.strand
            ) {safe_name}_m ON true
        """)
        
        select_cols.extend([
            f"{safe_name}_m.{safe_name}_start_matches",
            f"{safe_name}_m.{safe_name}_end_matches",
            f"{safe_name}_m.{safe_name}_both_matches"
        ])
    
    query = f"""
        SELECT 
            q.chr,
            q.start_pos,
            q.end_pos,
            q.strand,
            q.feature_name,
            q.score,
            {', '.join(select_cols)}
        FROM features q
        {' '.join(lateral_joins)}
        WHERE q.sample = '{sample_name}' 
          AND q.tool = '{query_tool}'
        ORDER BY q.chr, q.start_pos, q.end_pos
    """
    
    return con.execute(query).df()


def generate_ucsc_url(row, ucsc_session_url, flank=500):
    """Generate UCSC Genome Browser URL for a feature."""
    region = f"{row['chr']}:{max(0, row['start_pos']-flank)}-{row['end_pos']+flank}"
    return f"{ucsc_session_url}&position={region}"


def main():
    parser = argparse.ArgumentParser(
        description='Compare ncORF annotations across multiple tools',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example samplesheet format (TSV with header):
tool\tbed_path
ORFQuant\t/path/to/orfquant.bed
iRibo\t/path/to/iribo.bed
PRICE\t/path/to/price.bed
RiboTIE\t/path/to/ribotie.bed
        """
    )
    
    parser.add_argument(
        '-s', '--samplesheet',
        required=True,
        help='TSV file with columns: tool, bed_path'
    )
    
    parser.add_argument(
        '-n', '--sample-name',
        required=True,
        help='Sample name/identifier for this comparison'
    )
    
    parser.add_argument(
        '-o', '--outdir',
        default='consensus_results',
        help='Output directory (default: consensus_results)'
    )
    
    parser.add_argument(
        '-u', '--ucsc_session_url',
        default='hg38',
        help='UCSC Session URL (No position)'
    )
    
    parser.add_argument(
        '-f', '--flank',
        type=int,
        default=500,
        help='Flanking bases for UCSC links (default: 500)'
    )
    
    parser.add_argument(
        '--db',
        help='Path to DuckDB database (default: in-memory)'
    )
    
    args = parser.parse_args()
    
    # Setup output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Parse samplesheet
    print(f"Parsing samplesheet: {args.samplesheet}", file=sys.stderr)
    bed_files = parse_samplesheet(args.samplesheet)
    print(f"Found {len(bed_files)} tools: {', '.join(bed_files.keys())}", file=sys.stderr)
    
    # Initialize DuckDB
    db_path = args.db if args.db else ':memory:'
    con = duckdb.connect(db_path)
    
    try:
        # Ingest all BED files
        print(f"\nIngesting BED files for sample: {args.sample_name}", file=sys.stderr)
        ingest_beds(con, args.sample_name, bed_files)
        
        # Verify ingestion
        counts = con.execute(f"""
            SELECT tool, COUNT(*) as n_features
            FROM features
            WHERE sample = '{args.sample_name}'
            GROUP BY tool
        """).df()
        print("\nFeature counts per tool:", file=sys.stderr)
        print(counts.to_string(index=False), file=sys.stderr)
        
        # Compute matches for each tool
        print(f"\nComputing consensus matches...", file=sys.stderr)
        tools = list(bed_files.keys())
        
        for query_tool in tools:
            other_tools = [t for t in tools if t != query_tool]
            
            print(f"  Processing {query_tool}...", file=sys.stderr)
            df = compute_matches(con, args.sample_name, query_tool, other_tools)
            
            # Add UCSC links
            df['ucsc_url'] = df.apply(
                lambda row: generate_ucsc_url(row, args.ucsc_session_url, args.flank),
                axis=1
            )
            
            # Write output
            outfile = outdir / f"{args.sample_name}.{query_tool}.consensus.tsv"
            df.to_csv(outfile, sep='\t', index=False)
            print(f"    Wrote {len(df)} features to {outfile}", file=sys.stderr)
        
        # Generate summary statistics
        print(f"\nGenerating summary statistics...", file=sys.stderr)
        summary = []
        
        for query_tool in tools:
            df = pd.read_csv(outdir / f"{args.sample_name}.{query_tool}.consensus.tsv", sep='\t')
            
            stats = {
                'tool': query_tool,
                'n_features': len(df)
            }
            
            # Count features with matches in other tools
            for other_tool in [t for t in tools if t != query_tool]:
                safe_name = other_tool.replace('-', '_').replace('.', '_')
                
                stats[f'{other_tool}_start_any'] = (df[f'{safe_name}_start_matches'] > 0).sum()
                stats[f'{other_tool}_end_any'] = (df[f'{safe_name}_end_matches'] > 0).sum()
                stats[f'{other_tool}_both_any'] = (df[f'{safe_name}_both_matches'] > 0).sum()
            
            summary.append(stats)
        
        summary_df = pd.DataFrame(summary)
        summary_file = outdir / f"{args.sample_name}.summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"Wrote summary to {summary_file}", file=sys.stderr)
        
        print("\nDone!", file=sys.stderr)
        
    finally:
        con.close()


if __name__ == '__main__':
    main()