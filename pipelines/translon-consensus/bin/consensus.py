#!/usr/bin/env python3
"""
Consensus analysis for ncORF annotations across multiple tools.
Safer ingestion for GTF/BED and robust DuckDB CSV options.
"""

import argparse
import sys
from pathlib import Path
import gzip
import duckdb
import pandas as pd
from io import TextIOWrapper

# ---------- Utilities ----------

def is_gzip_file(path: Path) -> bool:
    """Check magic bytes to detect gzip (works even if extension is wrong)."""
    try:
        with open(path, 'rb') as fh:
            magic = fh.read(2)
            return magic == b'\x1f\x8b'
    except Exception:
        return False


def open_text_maybe_gz(path: Path):
    """Open file as text, transparently handling gzip if necessary."""
    if is_gzip_file(path):
        return TextIOWrapper(gzip.open(path, 'rb'), encoding='utf-8', errors='replace')
    else:
        return open(path, 'r', encoding='utf-8', errors='replace')


def find_first_data_line_with_tabs(path: Path, min_cols=9, max_lines=2000):
    """
    Return the first non-comment, non-empty line that contains at least `min_cols` tab-separated fields.
    Returns None if no such line in the first `max_lines`.
    """
    with open_text_maybe_gz(path) as fh:
        seen = 0
        for raw in fh:
            seen += 1
            if seen > max_lines:
                break
            line = raw.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                continue
            # Count fields by tabs
            n = line.count('\t') + 1
            if n >= min_cols:
                return line
        return None


# ---------- Samplesheet parsing ----------

def parse_samplesheet(samplesheet_path):
    """
    Parse samplesheet TSV with columns: tool, bed_path
    Returns: dict of {tool: bed_path}
    """
    df = pd.read_csv(samplesheet_path, sep='\t', header=0, dtype=str)
    
    if not {'tool', 'bed_path'}.issubset(df.columns):
        raise ValueError("Samplesheet must have columns: tool, bed_path")
    
    # Filter out missing files
    valid_entries = {}
    for _, row in df.iterrows():
        tool = str(row['tool'])
        bed_path = Path(row['bed_path'])
        
        if not bed_path.exists():
            print(f"Warning: BED file not found for {tool}: {bed_path}", file=sys.stderr)
            continue
        
        valid_entries[tool] = str(bed_path)
    
    if len(valid_entries) < 2:
        raise ValueError(f"Need at least 2 valid tool BED files, found {len(valid_entries)}")
    
    return valid_entries


# ---------- GTF ingestion (robust) ----------
from pathlib import Path
import sys


def find_first_data_line_with_tabs(path, min_cols=9, max_lines=2000):
    """
    Scan file for first non-comment line with >=min_cols tab-separated fields.
    Returns the line or None.
    """
    opener = open
    # Detect gzip magic manually (DuckDB can read gzip but our validator must too)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        import gzip
        opener = gzip.open

    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        count = 0
        for raw in fh:
            count += 1
            if count > max_lines:
                return None
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= min_cols:
                return line
    return None


def ingest_gtf_annotations(con, gtf_path):
    """
    Safely ingest a GENCODE-style GTF into DuckDB.
    Includes validation, compression detection, SQL brace escaping, and robust CSV options.
    """

    gtf_path = Path(gtf_path)
    print(f"Loading gencode annotations from: {gtf_path}", file=sys.stderr)

    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    # ---------- VALIDATION ----------
    data_line = find_first_data_line_with_tabs(gtf_path)
    if data_line is None:
        raise ValueError(
            f"No valid GTF data line found in {gtf_path}. "
            "Likely compressed, corrupt, or not a tab-separated GTF."
        )

    # ---------- TABLE CREATION ----------
    con.execute("""
        CREATE TABLE IF NOT EXISTS annotations (
            chr TEXT,
            source TEXT,
            feature_type TEXT,
            start_pos BIGINT,
            end_pos BIGINT,
            score TEXT,
            strand TEXT,
            frame TEXT,
            gene_id TEXT,
            gene_name TEXT,
            gene_biotype TEXT,
            transcript_id TEXT,
            transcript_name TEXT,
            transcript_biotype TEXT
        )
    """)

    # Escape path for SQL literal
    gtf_path_sql = str(gtf_path).replace("'", "''")

    # ---------- INSERT USING SAFE read_csv ----------
    #
    # IMPORTANT: columns={{ ... }} must use DOUBLE BRACES to escape Python f-string parsing.
    #
    con.execute(f"""
        INSERT INTO annotations
        SELECT
            seqname AS chr,
            "source",
            "feature" AS feature_type,
            "start" AS start_pos,
            "end" AS end_pos,
            "score",
            strand,
            "frame",
            regexp_extract(attributes, 'gene_id "([^"]+)"', 1) AS gene_id,
            regexp_extract(attributes, 'gene_name "([^"]+)"', 1) AS gene_name,
            regexp_extract(attributes, 'gene_type "([^"]+)"', 1) AS gene_biotype,
            regexp_extract(attributes, 'transcript_id "([^"]+)"', 1) AS transcript_id,
            regexp_extract(attributes, 'transcript_name "([^"]+)"', 1) AS transcript_name,
            regexp_extract(attributes, 'transcript_type "([^"]+)"', 1) AS transcript_biotype
        FROM read_csv(
            '{gtf_path_sql}',

            delim='\t',
            header=FALSE,
            quote='',              -- disable CSV quoting entirely
            escape='',             -- no escaping
            comment='#',           -- skip metadata/header

            ignore_errors=TRUE,    -- skip malformed records
            strict_mode=FALSE,     -- prevent dialect rejection
            null_padding=TRUE,     -- tolerate short rows
            sample_size=20000,     -- improve dialect inference

            columns={{             -- ESCAPED BRACES ({{ }}), not Python placeholders
                'seqname': 'VARCHAR',
                'source': 'VARCHAR',
                'feature': 'VARCHAR',
                'start': 'BIGINT',
                'end': 'BIGINT',
                'score': 'VARCHAR',
                'strand': 'VARCHAR',
                'frame': 'VARCHAR',
                'attributes': 'VARCHAR'
            }}
        )
        WHERE "feature" IN ('transcript', 'CDS', 'gene')
    """)

    print("GTF annotations loaded.", file=sys.stderr)

    # ---------- SUMMARY ----------
    counts = con.execute("""
        SELECT feature_type, COUNT(*) AS n
        FROM annotations
        GROUP BY feature_type
        ORDER BY n DESC
    """).df()

    print(counts.to_string(index=False), file=sys.stderr)

    # ---------- INDEXES ----------
    con.execute("CREATE INDEX IF NOT EXISTS idx_ann_loc ON annotations(chr, strand, start_pos, end_pos)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_ann_tx  ON annotations(transcript_id)")

# ---------- BED ingestion (robust) ----------

def ingest_beds(con, sample_name, bed_files_dict):
    """
    Load all BED files into DuckDB features table, consolidating duplicate features.

    Uses read_csv_auto with tolerant options; handles gzipped BEDs and deduplicates
    using ROW_NUMBER() over partition (DuckDB-friendly).
    """
    # Create features table
    con.execute('''
        CREATE TABLE IF NOT EXISTS features (
            sample TEXT,
            tool TEXT,
            chr TEXT,
            chromStart BIGINT,
            chromEnd BIGINT,
            feature_name TEXT,
            score REAL,
            strand TEXT,
            start_pos BIGINT,
            end_pos BIGINT,
            itemRgb TEXT,
            blockCount INTEGER,
            blockSizes TEXT,
            blockStarts TEXT
        )
    ''')

    for tool, bed_path in bed_files_dict.items():
        bed_path = Path(bed_path)
        print(f"Ingesting {tool}: {bed_path}", file=sys.stderr)

        if not bed_path.exists():
            print(f"  ERROR: {bed_path} not found, skipping {tool}", file=sys.stderr)
            continue

        # Read first non-empty, non-comment line to count fields.
        # Support gz as well.
        first_data = None
        with open_text_maybe_gz(bed_path) as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                    continue
                first_data = line
                break

        if first_data is None:
            print(f"  WARNING: Could not read a data line from {bed_path}, skipping.", file=sys.stderr)
            continue

        n_fields = len(first_data.split('\t'))
        if n_fields >= 12:
            columns = [
                ('chr', 'VARCHAR'),
                ('chromStart', 'BIGINT'),  # Keep original for reference
                ('chromEnd', 'BIGINT'),    # Keep original for reference
                ('feature_name', 'VARCHAR'),
                ('score', 'REAL'),
                ('strand', 'VARCHAR'),
                ('start_pos', 'BIGINT'),   # Use thickStart as actual start
                ('end_pos', 'BIGINT'),     # Use thickEnd as actual end
                ('itemRgb', 'VARCHAR'),
                ('blockCount', 'INTEGER'),
                ('blockSizes', 'VARCHAR'),
                ('blockStarts', 'VARCHAR')
            ]
        elif n_fields >= 6:
            # BED6: chr, start, end, name, score, strand
            # For BED6, chromStart/chromEnd are same as start_pos/end_pos (no thick/thin distinction)
            columns = [
                ('chr', 'VARCHAR'),
                ('chromStart', 'BIGINT'),   # Column 2 (same as start_pos for BED6)
                ('chromEnd', 'BIGINT'),     # Column 3 (same as end_pos for BED6)
                ('feature_name', 'VARCHAR'),
                ('score', 'REAL'),
                ('strand', 'VARCHAR')
            ]
        elif n_fields >= 3:
            # BED3: chr, start, end
            # For BED3, chromStart/chromEnd are same as start_pos/end_pos
            columns = [
                ('chr', 'VARCHAR'),
                ('chromStart', 'BIGINT'),   # Column 2 (same as start_pos for BED3)
                ('chromEnd', 'BIGINT')      # Column 3 (same as end_pos for BED3)
            ]
        else:
            raise ValueError(f"Invalid BED format in {bed_path}: only {n_fields} fields")

        # Build column spec for read_csv_auto
        col_spec = {f'col{i}': dtype for i, (_, dtype) in enumerate(columns.items())} if isinstance(columns, dict) else {f'col{i}': dtype for i, (_, dtype) in enumerate(columns)}
        # But above we used list of tuples; handle accordingly
        if isinstance(columns, list):
            col_spec = {f'col{i}': dtype for i, (_, dtype) in enumerate(columns)}
            col_names = [name for name, _ in columns]
        else:
            col_names = list(columns.keys())

        # Prepare select columns mapping read_csv_auto colN -> final names
        select_cols = ', '.join([f'col{i} as {name}' for i, name in enumerate(col_names)])

        # For BED6/BED3, add start_pos/end_pos as aliases for chromStart/chromEnd
        extra_cols = []
        if 'chromStart' in col_names and 'start_pos' not in col_names:
            extra_cols.append('chromStart as start_pos')
        if 'chromEnd' in col_names and 'end_pos' not in col_names:
            extra_cols.append('chromEnd as end_pos')

        all_select_cols = select_cols
        if extra_cols:
            all_select_cols = select_cols + ', ' + ', '.join(extra_cols)

        bed_path_sql = str(bed_path).replace("'", "''")
        sample_sql = str(sample_name).replace("'", "''")
        tool_sql = str(tool).replace("'", "''")

        try:
            # Create temp table
            con.execute(f"""
                CREATE OR REPLACE TEMP TABLE temp_features AS
                SELECT
                    '{sample_sql}' as sample,
                    '{tool_sql}' as tool,
                    {all_select_cols}
                FROM read_csv_auto('{bed_path_sql}',
                    delim='\\t',
                    header=false,
                    ignore_errors=true,
                    null_padding=true,
                    sample_size=20000,
                    columns={{{', '.join([f"'{k}': '{v}'" for k, v in col_spec.items()])}}}
                )
            """)

            # Count before dedup
            count_before = con.execute("SELECT COUNT(*) FROM temp_features").fetchone()[0]

            # Deduplicate using row_number() partition (DuckDB supports window functions)
            # Explicitly specify column order to match table schema
            table_col_order = [
                'chr', 'chromStart', 'chromEnd', 'feature_name', 'score',
                'strand', 'start_pos', 'end_pos', 'itemRgb', 'blockCount',
                'blockSizes', 'blockStarts'
            ]
            # Build list of available columns (includes both col_names and extra_cols)
            available_cols = set(col_names)
            if 'chromStart' in col_names and 'start_pos' not in col_names:
                available_cols.add('start_pos')
            if 'chromEnd' in col_names and 'end_pos' not in col_names:
                available_cols.add('end_pos')

            # Only include columns that exist in temp_features table
            insert_cols = ', '.join([c for c in table_col_order if c in available_cols])
            con.execute(f"""
                INSERT INTO features (sample, tool, {insert_cols})
                SELECT sample, tool, {insert_cols} FROM (
                    SELECT *,
                        ROW_NUMBER() OVER (
                            PARTITION BY chr, start_pos, end_pos, COALESCE(strand,''), COALESCE(feature_name,'')
                            ORDER BY chr, start_pos, end_pos
                        ) as rn
                    FROM temp_features
                ) t
                WHERE rn = 1
            """)

            count_after = con.execute(f"""
                SELECT COUNT(*) FROM features
                WHERE sample = '{sample_sql}' AND tool = '{tool_sql}'
            """).fetchone()[0]

            duplicates_removed = count_before - count_after
            if duplicates_removed > 0:
                print(f"  → Removed {duplicates_removed} duplicate features ({count_before} → {count_after})", file=sys.stderr)
            else:
                print(f"  → Inserted {count_after} features for {tool}", file=sys.stderr)

        except Exception as e:
            print(f"Error ingesting {tool} from {bed_path}: {e}", file=sys.stderr)
            # ensure temp table is dropped if exists
            try:
                con.execute("DROP TABLE IF EXISTS temp_features")
            except Exception:
                pass
            raise
        finally:
            # Try to clean up temp table
            try:
                con.execute("DROP TABLE IF EXISTS temp_features")
            except Exception:
                pass


# ---------- Matching / Annotation / Output ----------

def compute_matches(con, sample_name, query_tool, other_tools):
    """
    For each feature in query_tool, count exact matches and overlaps in all other tools.
    Uses separate queries for start/end exact matches, plus interval overlaps.
    """
    lateral_joins = []
    select_cols = []

    for other_tool in other_tools:
        safe_name = other_tool.replace('-', '_').replace('.', '_')

        # Exact start matches
        lateral_joins.append(f"""
            LEFT JOIN LATERAL (
                SELECT
                    COALESCE(SUM(CASE WHEN start_pos = q.start_pos THEN 1 ELSE 0 END), 0) as {safe_name}_start_matches
                FROM features
                WHERE sample = '{sample_name.replace("'", "''")}'
                  AND tool = '{other_tool.replace("'", "''")}'
                  AND chr = q.chr
                  AND COALESCE(strand,'') = COALESCE(q.strand,'')
            ) {safe_name}_start ON true
        """)

        # Exact end matches
        lateral_joins.append(f"""
            LEFT JOIN LATERAL (
                SELECT
                    COALESCE(SUM(CASE WHEN end_pos = q.end_pos THEN 1 ELSE 0 END), 0) as {safe_name}_end_matches
                FROM features
                WHERE sample = '{sample_name.replace("'", "''")}'
                  AND tool = '{other_tool.replace("'", "''")}'
                  AND chr = q.chr
                  AND COALESCE(strand,'') = COALESCE(q.strand,'')
            ) {safe_name}_end ON true
        """)

        # Interval overlaps (any overlap, not just exact)
        lateral_joins.append(f"""
            LEFT JOIN LATERAL (
                SELECT
                    COALESCE(COUNT(*), 0) as {safe_name}_overlap_count
                FROM features
                WHERE sample = '{sample_name.replace("'", "''")}'
                  AND tool = '{other_tool.replace("'", "''")}'
                  AND chr = q.chr
                  AND COALESCE(strand,'') = COALESCE(q.strand,'')
                  AND start_pos < q.end_pos
                  AND end_pos > q.start_pos
            ) {safe_name}_overlap ON true
        """)

        select_cols.extend([
            f"{safe_name}_start.{safe_name}_start_matches",
            f"{safe_name}_end.{safe_name}_end_matches",
            # Compute both as the minimum of start and end matches
            f"LEAST({safe_name}_start.{safe_name}_start_matches, {safe_name}_end.{safe_name}_end_matches) as {safe_name}_both_matches",
            f"{safe_name}_overlap.{safe_name}_overlap_count"
        ])

    query = f"""
        SELECT
            q.chr,
            q.start_pos,
            q.end_pos,
            q.strand,
            q.feature_name,
            q.score,
            q.blockCount,
            q.blockSizes,
            -- Calculate CDS length by summing blockSizes (exon lengths)
            -- If blockSizes is NULL or empty, fall back to end_pos - start_pos
            CASE
                WHEN q.blockSizes IS NOT NULL AND q.blockSizes != '' THEN
                    (SELECT SUM(CAST(unnest(string_split(q.blockSizes, ',')) AS INTEGER)))
                ELSE
                    q.end_pos - q.start_pos
            END AS cds_length,
            {', '.join(select_cols)}
        FROM features q
        {' '.join(lateral_joins)}
        WHERE q.sample = '{sample_name.replace("'", "''")}'
          AND q.tool = '{query_tool.replace("'", "''")}'
        ORDER BY q.chr, q.start_pos, q.end_pos
    """

    return con.execute(query).df()


def annotate_features(con, sample_name, query_tool):
    """
    Annotate features with gencode transcript biotype and CDS position context.
    """
    query = f"""
    WITH feature_annotations AS (
        SELECT
            f.*,
            t.transcript_id,
            t.transcript_biotype,
            t.gene_name,
            CASE
                WHEN t.transcript_id IS NULL THEN 'intergenic'
                ELSE
                    CASE
                        WHEN EXISTS (
                            SELECT 1 FROM annotations cds
                            WHERE cds.transcript_id = t.transcript_id
                            AND cds.feature_type = 'CDS'
                        ) THEN
                            CASE
                                WHEN f.end_pos < (
                                    SELECT MIN(start_pos) FROM annotations cds
                                    WHERE cds.transcript_id = t.transcript_id
                                    AND cds.feature_type = 'CDS'
                                ) THEN '5prime_UTR'
                                WHEN f.start_pos > (
                                    SELECT MAX(end_pos) FROM annotations cds
                                    WHERE cds.transcript_id = t.transcript_id
                                    AND cds.feature_type = 'CDS'
                                ) THEN '3prime_UTR'
                                WHEN f.start_pos < (
                                    SELECT MIN(start_pos) FROM annotations cds
                                    WHERE cds.transcript_id = t.transcript_id
                                    AND cds.feature_type = 'CDS'
                                ) AND f.end_pos > (
                                    SELECT MIN(start_pos) FROM annotations cds
                                    WHERE cds.transcript_id = t.transcript_id
                                    AND cds.feature_type = 'CDS'
                                ) THEN 'overlaps_CDS_5prime'
                                WHEN f.start_pos < (
                                    SELECT MAX(end_pos) FROM annotations cds
                                    WHERE cds.transcript_id = t.transcript_id
                                    AND cds.feature_type = 'CDS'
                                ) AND f.end_pos > (
                                    SELECT MAX(end_pos) FROM annotations cds
                                    WHERE cds.transcript_id = t.transcript_id
                                    AND cds.feature_type = 'CDS'
                                ) THEN 'overlaps_CDS_3prime'
                                ELSE 'within_CDS'
                            END
                        ELSE 'within_transcript_no_CDS'
                    END
            END as cds_context
        FROM features f
        LEFT JOIN annotations t ON
            f.chr = t.chr
            AND COALESCE(f.strand,'') = COALESCE(t.strand,'')
            AND t.feature_type = 'transcript'
            AND NOT (f.end_pos < t.start_pos OR f.start_pos > t.end_pos)
        WHERE f.sample = '{sample_name.replace("'", "''")}'
        AND f.tool = '{query_tool.replace("'", "''")}'
    )
    SELECT
        chr,
        start_pos,
        end_pos,
        strand,
        feature_name,
        score,
        COALESCE(transcript_id, 'none') as transcript_id,
        COALESCE(transcript_biotype, 'intergenic') as rna_biotype,
        COALESCE(gene_name, 'intergenic') as gene_name,
        cds_context
    FROM feature_annotations
    """

    return con.execute(query).df()


def generate_ucsc_url(row, ucsc_session_url, flank=500):
    """Generate UCSC Genome Browser URL for a feature."""
    region = f"{row['chr']}:{max(0, int(row['start_pos'])-int(flank))}-{int(row['end_pos'])+int(flank)}"
    return f"{ucsc_session_url}&position={region}"


# ---------- Main ----------

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
    
    parser.add_argument('-s', '--samplesheet', required=True, help='TSV file with columns: tool, bed_path')
    parser.add_argument('-n', '--sample-name', required=True, help='Sample name/identifier for this comparison')
    parser.add_argument('-o', '--outdir', default='consensus_results', help='Output directory (default: consensus_results)')
    parser.add_argument('-u', '--ucsc_session_url', default='hg38', help='UCSC Session URL (No position)')
    parser.add_argument('-f', '--flank', type=int, default=500, help='Flanking bases for UCSC links (default: 500)')
    parser.add_argument('-g', '--gencode-gtf', help='Path to gencode GTF file for annotation (can be .gz). If not provided, features will not be annotated with biotype/CDS context.')
    parser.add_argument('--db', help='Path to DuckDB database (default: in-memory)')
    
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
        # Ingest BED files
        print(f"\nIngesting BED files for sample: {args.sample_name}", file=sys.stderr)
        ingest_beds(con, args.sample_name, bed_files)

        # Load GTF annotations if provided
        has_annotations = False
        if args.gencode_gtf:
            gtf_path = Path(args.gencode_gtf)
            if gtf_path.exists():
                print(f"\nLoading gencode annotations...", file=sys.stderr)
                ingest_gtf_annotations(con, str(gtf_path))
                has_annotations = True
            else:
                print(f"Warning: GTF file not found: {gtf_path}", file=sys.stderr)
                print("  Proceeding without annotation...", file=sys.stderr)

        # Verify ingestion
        counts = con.execute(f"""
            SELECT tool, COUNT(*) as n_features
            FROM features
            WHERE sample = '{args.sample_name.replace("'", "''")}'
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

            # Add annotations if available
            if has_annotations:
                print(f"    Adding gencode annotations...", file=sys.stderr)
                annot_df = annotate_features(con, args.sample_name, query_tool)

                # Merge annotations into consensus results
                df = df.merge(
                    annot_df[['chr', 'start_pos', 'end_pos', 'strand',
                             'transcript_id', 'rna_biotype', 'gene_name', 'cds_context']],
                    on=['chr', 'start_pos', 'end_pos', 'strand'],
                    how='left'
                )

                # Fill missing values for features without annotation
                df['transcript_id'] = df['transcript_id'].fillna('none')
                df['rna_biotype'] = df['rna_biotype'].fillna('intergenic')
                df['gene_name'] = df['gene_name'].fillna('intergenic')
                df['cds_context'] = df['cds_context'].fillna('intergenic')

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
                
                stats[f'{other_tool}_start_any'] = int((df[f'{safe_name}_start_matches'] > 0).sum())
                stats[f'{other_tool}_end_any'] = int((df[f'{safe_name}_end_matches'] > 0).sum())
                stats[f'{other_tool}_both_any'] = int((df[f'{safe_name}_both_matches'] > 0).sum())
            
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
