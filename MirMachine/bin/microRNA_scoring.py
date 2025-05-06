#!/usr/bin/env python
import pandas as pd
import argparse
import re

def convert_heatmap_to_flat_tsv(input_path, output_path):
    """Convert the alternating format heatmap to a flat TSV with one entry per line"""
    flat_data = []
    current_species = None
    
    with open(input_path, 'r') as infile:
        # Read header
        header = next(infile).strip().split(',')
        
        for line in infile:
            line = line.strip()
            if not line:
                continue
            
            # Species lines have few or no commas
            if line.count(',') <= 1:
                current_species = line
                continue
            
            # Data lines have multiple fields
            parts = line.split(',')
            
            # Create a flat entry with species as first column
            row = [current_species] + parts
            flat_data.append(row)
    
    # Write the flat TSV
    with open(output_path, 'w') as outfile:
        # Add species to the header
        outfile.write("full_species\t" + "\t".join(header) + "\n")
        
        # Write all rows
        for row in flat_data:
            outfile.write("\t".join(row) + "\n")
    
    print(f"Converted {len(flat_data)} entries to flat TSV format: {output_path}")
    return len(flat_data)

def parse_metadata(metadata_path):
    """Extract total families searched for each species from metadata file"""
    species_to_families = {}
    
    with open(metadata_path, 'r') as f:
        for line in f:
            # Extract total families searched and species
            match = re.search(r'#,Total,families,searched:,(\d+),#,Species:,(.+)', line.strip())
            if match:
                families = int(match.group(1))
                species = match.group(2)
                species_to_families[species] = families
    print(f"Number of species in metadata: {len(species_to_families)}")
    return species_to_families

def parse_flat_tsv(flat_tsv_path, metadata_dict):
    """Parse the flattened TSV file and calculate microRNA scores"""
    # Read the flat TSV
    df = pd.read_csv(flat_tsv_path, sep='\t')
    
    # Extract species and assembly_accession
    def extract_base_species(full_species):
        if '_GCA_' in full_species:
            return full_species.split('_GCA_')[0]
        elif '_GCF_' in full_species:
            return full_species.split('_GCF_')[0]
        return full_species
    
    def extract_accession(full_species):
        if '_GCA_' in full_species:
            return 'GCA_' + full_species.split('_GCA_')[1]
        elif '_GCF_' in full_species:
            return 'GCF_' + full_species.split('_GCF_')[1]
        return None
    
    # Add columns for base species and accession
    df['species'] = df['full_species'].apply(extract_base_species)
    df['assembly_accession'] = df['full_species'].apply(extract_accession)
    
    # Print statistics
    print(f"Read {len(df)} entries from flat TSV")
    print(f"Found {len(df['full_species'].unique())} unique full species")
    print(f"Found {len(df['species'].unique())} unique base species")
    print(f"Found {len(df['assembly_accession'].dropna().unique())} unique assembly accessions")
    
    # Process the data by species and assembly
    results = {}
    
    # Group by species and assembly accession
    for (species, accession), group in df.groupby(['species', 'assembly_accession']):
        # Skip entries with no accession
        if pd.isna(accession):
            continue
            
        # Get total families searched for this species
        total_families = metadata_dict.get(species, 0)
        if total_families == 0:
            # Use the number of families in the data as a fallback
            unique_families = group['family'].dropna().unique()
            total_families = len(unique_families)
            if total_families == 0:
                total_families = 1  # Avoid division by zero
            print(f"Warning: No total families found for {species}, using {total_families} from data")
            
        # Get mode/analysis node
        analysis_node = group['mode'].iloc[0] if 'mode' in group.columns and not group['mode'].empty else "Unknown"
        
        # Filtered data (has filtered value)
        filtered_df = group.dropna(subset=['filtered'])
        filtered_count = len(filtered_df)
        filtered_score = (filtered_count / total_families) * 100
        
        # Unfiltered data (has tgff value)
        unfiltered_df = group.dropna(subset=['tgff'])
        unfiltered_count = len(unfiltered_df)
        unfiltered_score = (unfiltered_count / total_families) * 100
        
        # Families that were filtered out (has tgff but no filtered)
        filtered_out = group[group['filtered'].isna() & group['tgff'].notna()]
        filtered_out_families = ", ".join(filtered_out['family'].dropna())
        
        # Families with no hits (no tgff and no filtered)
        no_hits = group[group['tgff'].isna() & group['filtered'].isna()]
        filtered_no_hits = ", ".join(no_hits['family'].dropna())
        
        # Store results for this assembly
        results[f"{species}_{accession}"] = {
            'species': species,
            'assembly_accession': accession,
            'analysis_node': analysis_node,
            'filtered_total_count': filtered_count,
            'filtered_microRNA_score': round(filtered_score, 2),
            'filtered_no_hits': filtered_no_hits,
            'unfiltered_total_count': unfiltered_count,
            'unfiltered_microRNA_score': round(unfiltered_score, 2),
            'filtered_out_families': filtered_out_families
        }
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(list(results.values()))
    
    # Reorder columns
    column_order = ["species", "assembly_accession", "analysis_node", "filtered_total_count", 
                   "filtered_microRNA_score", "filtered_no_hits", 
                   "unfiltered_total_count", "unfiltered_microRNA_score",
                   "filtered_out_families"]
    
    results_df = results_df.reindex(columns=column_order)
    
    print(f"Successfully processed {len(results_df)} assembly accessions")
    return results_df

def main():
    parser = argparse.ArgumentParser(description='Score microRNA sequences')
    parser.add_argument('-i', '--input', required=True, help='Path to the miRNA heatmap CSV file')
    parser.add_argument('-m', '--metadata', required=True, help='Path to the miRNA metadata file')
    parser.add_argument('-o', '--output', required=True, help='Path to output results CSV')
    parser.add_argument('--json', action='store_true', help='Also output results as JSON')
    args = parser.parse_args()
    
    # Parse input files
    metadata_dict = parse_metadata(args.metadata)

    # Convert heatmap to flat CSV if needed
    if args.input.endswith('.csv'):
        flat_tsv_path = args.input.replace('.csv', '_flat.tsv')
        convert_heatmap_to_flat_tsv(args.input, flat_tsv_path)
        args.input = flat_tsv_path
    
    results_df = parse_flat_tsv(args.input, metadata_dict)

    print(results_df.head())
    # Write output
    results_df.to_csv(args.output, index=False, sep='\t' if args.output.endswith('.tsv') else ',')
    
    if args.json:
        json_path = args.output.replace('.csv', '.json').replace('.tsv', '.json')
        results_df.to_json(json_path, orient='records', indent=4)
    
    print(f"Successfully processed {len(results_df)} species")

if __name__ == '__main__':
    main()