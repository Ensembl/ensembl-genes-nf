#!/usr/bin/env python



import pandas as pd
import re
import sys
import argparse

def process_mirna_file(input_file, output_file):
    """Process a single microRNA heatmap file and create a standardized TSV output."""
    
    # Read the file content
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Extract metadata
    species_match = re.search(r'# Species: (.+)', content)
    species = species_match.group(1) if species_match else "Unknown"
    
    # Extract assembly accession from the genome file line
    genome_match = re.search(r'# Genome file: (.+)\.fa', content)
    assembly_accession = genome_match.group(1) if genome_match else "Unknown"
    
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
        if line.startswith('species,query_node,family,node,total_hits,filtered_hits'):
            data_start = i
            break
    
    if data_start == 0:
        print("Error: Could not find data header in the file.")
        return
    
    # Create a dataframe from the CSV data
    data_lines = lines[data_start:]
    data_text = '\n'.join(data_lines)
    df = pd.read_csv(pd.StringIO(data_text))
    
    # Calculate filtered and unfiltered counts
    filtered_df = df[df['filtered_hits'].notna()]
    filtered_count = len(filtered_df)
    filtered_score = (filtered_count / total_families) * 100
    
    unfiltered_df = df[df['total_hits'].notna()]
    unfiltered_count = len(unfiltered_df)
    unfiltered_score = (unfiltered_count / total_families) * 100
    
    # Find families with no hits
    no_hits_df = df[(df['total_hits'].isna()) & (df['filtered_hits'].isna())]
    filtered_no_hits = ", ".join(no_hits_df['family'].tolist())
    
    # Find families that were filtered out (have total hits but no filtered hits)
    filtered_out_df = df[df['filtered_hits'].isna() & df['total_hits'].notna()]
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
    
    # Convert to DataFrame and save as TSV
    result_df = pd.DataFrame([result])
    
    # Reorder columns
    column_order = ["species", "assembly_accession", "analysis_node", "filtered_total_count", 
                   "filtered_microRNA_score", "filtered_no_hits", 
                   "unfiltered_total_count", "unfiltered_microRNA_score",
                   "filtered_out_families"]
    
    result_df = result_df.reindex(columns=column_order)
    
    # Write to TSV
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed {species} ({assembly_accession}) and saved to {output_file}")
    print(f"Filtered score: {filtered_score:.2f}%, Unfiltered score: {unfiltered_score:.2f}%")

def main():
    parser = argparse.ArgumentParser(description='Process microRNA data for a single assembly')
    parser.add_argument('-i', '--input', required=True, help='Input microRNA heatmap file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()
    
    process_mirna_file(args.input, args.output)

if __name__ == "__main__":
    main()