'''
Scoring Methodology

Filtered microRNA Score: Percentage of microRNA families passing quality filters (filtered hits/total families × 100)
Unfiltered microRNA Score: Percentage of families with any matches before filtering (unfiltered hits/total families × 100)

Key Metrics

Species: Organism analyzed
Assembly Accession: Genome assembly ID (GCA/GCF format)
Analysis Node: Computational node used
Filtered Total Count: Number of families passing filters
Filtered No Hits: Families with no matches
Unfiltered Total Count: Families with any matches
Filtered Out Families: Families detected but filtered out

'''


import pandas as pd
import re
import os
import sys
import glob
from io import StringIO  # Add this import

def process_mirna_file(input_file):
    """Process a single microRNA heatmap file and return a standardized row."""
    
    # Read the file content
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Extract metadata
    species_match = re.search(r'# Species: (.+)', content)
    species = species_match.group(1) if species_match else "Unknown"
    
    # Extract assembly accession from the genome file line
    genome_match = re.search(r'# Genome file: (.+)\.fa', content)
    assembly_accession = None
    if genome_match:
        accession_str = genome_match.group(1)
        # Check if it's in GCA/GCF format
        if "GCA_" in accession_str or "GCF_" in accession_str:
            assembly_accession = accession_str
        else:
            # Try to extract from filename if not in genome line
            filename = os.path.basename(input_file)
            acc_match = re.search(r'(GCA_\d+\.\d+|GCF_\d+\.\d+)', filename)
            if acc_match:
                assembly_accession = acc_match.group(1)
            else:
                assembly_accession = "Unknown"
    else:
        # Try to extract from filename if not in genome line
        filename = os.path.basename(input_file)
        acc_match = re.search(r'(GCA_\d+\.\d+|GCF_\d+\.\d+)', filename)
        if acc_match:
            assembly_accession = acc_match.group(1)
        else:
            assembly_accession = "Unknown"
    
    # Get the analysis node
    node_match = re.search(r'# Node: (.+)', content)
    analysis_node = node_match.group(1) if node_match else "Unknown"
    
    # Get the total families searched
    families_match = re.search(r'# Total families searched: (\d+)', content)
    total_families = int(families_match.group(1)) if families_match else 0
    
    # Parse the data portion into a DataFrame
    lines = content.split('\n')
    data_start = 0
    for i, line in enumerate(lines):
        if 'species' in line and 'family' in line and ('total_hits' in line or 'tgff' in line):
            data_start = i
            break
    
    if data_start == 0:
        print(f"Error: Could not find data header in the file {input_file}", file=sys.stderr)
        return None
    
    # Create a dataframe from the CSV data
    data_lines = lines[data_start:]
    data_text = '\n'.join(data_lines)
    
    try:
        df = pd.read_csv(StringIO(data_text))  # Fixed here
    except Exception as e:
        print(f"Error parsing CSV data from {input_file}: {e}", file=sys.stderr)
        return None
    
    # Handle column name differences (some files use tgff instead of total_hits)
    if 'total_hits' in df.columns and 'filtered_hits' in df.columns:
        total_col = 'total_hits'
        filtered_col = 'filtered_hits'
    elif 'tgff' in df.columns and 'filtered' in df.columns:
        total_col = 'tgff'
        filtered_col = 'filtered'
    else:
        print(f"Error: Could not determine column names in {input_file}", file=sys.stderr)
        return None
    
    # Calculate filtered and unfiltered counts
    filtered_df = df[df[filtered_col].notna()]
    filtered_count = len(filtered_df)
    filtered_score = (filtered_count / total_families) * 100 if total_families > 0 else 0
    
    unfiltered_df = df[df[total_col].notna()]
    unfiltered_count = len(unfiltered_df)
    unfiltered_score = (unfiltered_count / total_families) * 100 if total_families > 0 else 0
    
    # Find families with no hits
    no_hits_df = df[(df[total_col].isna()) & (df[filtered_col].isna())]
    filtered_no_hits = ", ".join(no_hits_df['family'].tolist())
    
    # Find families that were filtered out (have total hits but no filtered hits)
    filtered_out_df = df[df[filtered_col].isna() & df[total_col].notna()]
    filtered_out_families = ", ".join(filtered_out_df['family'].tolist())
    
    # Create the result dictionary
    result = {
        'species': species,
        'assembly_accession': assembly_accession,
        'analysis_node': analysis_node,
        'filtered_total_count': filtered_count,
        'filtered_microRNA_score': round(filtered_score, 2),
        'filtered_no_hits': filtered_no_hits,
        'unfiltered_total_count': unfiltered_count,
        'unfiltered_microRNA_score': round(unfiltered_score, 2),
        'filtered_out_families': filtered_out_families
    }
    
    return result

# Process all files from command line arguments
def process_files(file_pattern):
    results = []
    files = glob.glob(file_pattern)
    total_files = len(files)
    
    for i, file in enumerate(files):
        if (i+1) % 100 == 0:
            print(f"Processing file {i+1}/{total_files}: {file}", file=sys.stderr)
        
        result = process_mirna_file(file)
        if result is not None:
            results.append(result)
    
    # Convert to DataFrame
    if results:
        result_df = pd.DataFrame(results)
        
        # Reorder columns
        column_order = ["species", "assembly_accession", "analysis_node", "filtered_total_count", 
                       "filtered_microRNA_score", "filtered_no_hits", 
                       "unfiltered_total_count", "unfiltered_microRNA_score",
                       "filtered_out_families"]
        
        result_df = result_df.reindex(columns=column_order)
        
        # Write to stdout as TSV
        result_df.to_csv(sys.stdout, sep='\t', index=False)
        print(f"Processed {len(results)} files successfully", file=sys.stderr)
    else:
        print("No results were produced", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python process_mirna.py <file_pattern>", file=sys.stderr)
        sys.exit(1)
    
    process_files(sys.argv[1])
