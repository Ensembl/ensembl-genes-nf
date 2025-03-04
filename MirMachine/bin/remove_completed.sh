#!/bin/bash

# Usage: ./script.sh tsv_file directory_to_search
# Example: ./script.sh genomes.tsv /path/to/genomes

if [ $# -ne 2 ]; then
    echo "Usage: $0 <tsv_file> <directory_to_search>"
    exit 1
fi

tsv_file="$1"
search_dir="$2"

# Check if files exist
if [ ! -f "$tsv_file" ]; then
    echo "Error: TSV file '$tsv_file' not found"
    exit 1
fi

if [ ! -d "$search_dir" ]; then
    echo "Error: Directory '$search_dir' not found"
    exit 1
fi

# Create a temporary file to store filtered accessions
temp_file=$(mktemp)

# Process the TSV file
while IFS=$'\t' read -r accession rest_of_line; do
    # Skip header line if needed (uncomment if your TSV has a header)
    # [ "$accession" = "GCA_Accession" ] && continue
    
    # Check if a directory with this accession name exists
    if [ ! -d "$search_dir/$accession" ]; then
        # If no directory exists, add this line to our filtered list
        echo -e "$accession\t$rest_of_line" >> "$temp_file"
    fi
done < "$tsv_file"

# Report results
total_genomes=$(wc -l < "$tsv_file")
filtered_genomes=$(wc -l < "$temp_file")

echo "Total genomes in TSV: $total_genomes"
echo "Genomes without directories: $filtered_genomes"
echo "Filtered results saved to: $temp_file"

# Optionally, you can display the first few filtered accessions
echo "First 5 filtered accessions:"
head -n 5 "$temp_file"