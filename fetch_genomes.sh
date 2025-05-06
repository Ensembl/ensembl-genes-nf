#!/bin/bash

# Script to iterate over a genome list file and run rapid_fetch.py for each entry
# Optimized for CSV format with species names that may contain multiple words
#
# Usage:
#   ./fetch_multiple_genomes.sh -f <genome_list_file> -o <output_directory> [-g <gca_column>] [-s <species_columns>]
#
# Options:
#   -f <genome_list_file>   File containing genome information in CSV format
#   -o <output_directory>   Directory to save fetched genomes
#   -g <gca_column>         Column number (1-based) containing the GCA ID (default: auto-detect)
#   -s <species_columns>    Number of columns to use for species name (default: auto-detect)

# Default values
OUTPUT_DIR="./genomes"
GENOME_FILE=""
GCA_COLUMN=0  # 0 means auto-detect
SPECIES_COLUMNS=0  # 0 means auto-detect

# Parse command line arguments
while getopts "f:o:g:s:" opt; do
  case ${opt} in
    f )
      GENOME_FILE=$OPTARG
      ;;
    o )
      OUTPUT_DIR=$OPTARG
      ;;
    g )
      GCA_COLUMN=$OPTARG
      ;;
    s )
      SPECIES_COLUMNS=$OPTARG
      ;;
    \? )
      echo "Usage: $0 -f <genome_list_file> -o <output_directory> [-g <gca_column>] [-s <species_columns>]"
      exit 1
      ;;
  esac
done

# Check if genome file was provided
if [ -z "$GENOME_FILE" ]; then
  echo "Error: Genome list file (-f) is required"
  echo "Usage: $0 -f <genome_list_file> -o <output_directory> [-g <gca_column>] [-s <species_columns>]"
  exit 1
fi

# Check if the genome file exists
if [ ! -f "$GENOME_FILE" ]; then
  echo "Error: Genome list file '$GENOME_FILE' does not exist"
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Initialize counters
TOTAL=0
SUCCESSFUL=0
FAILED=0

# Process each line in the genome file
while IFS= read -r line; do
  # Skip empty lines
  if [ -z "$line" ]; then
    continue
  fi
  
  # Remove any quotes and split the line by comma
  line=$(echo "$line" | tr -d '"')
  IFS=',' read -ra FIELDS <<< "$line"
  
  # Skip lines with too few fields
  if [ ${#FIELDS[@]} -lt 2 ]; then
    echo "Warning: Skipping line with insufficient fields: $line"
    continue
  fi
  
  # Auto-detect GCA column if not specified
  if [ $GCA_COLUMN -eq 0 ]; then
    for i in "${!FIELDS[@]}"; do
      if [[ ${FIELDS[$i]} == GCA_* ]]; then
        GCA_COLUMN=$((i + 1))  # Convert to 1-based indexing
        break
      fi
    done
    
    if [ $GCA_COLUMN -eq 0 ]; then
      echo "Warning: No GCA ID found in line: $line"
      continue
    fi
  fi
  
  # Get GCA ID from the specified column (convert from 1-based to 0-based)
  GCA=${FIELDS[$((GCA_COLUMN - 1))]}
  
  # Skip if no valid GCA ID
  if [[ ! $GCA =~ ^GCA_ ]]; then
    echo "Warning: Invalid GCA ID format in column $GCA_COLUMN: $GCA"
    continue
  fi
  
  # Auto-detect species columns if not specified
  # Strategy: Everything before the GCA column is considered part of the species name
  if [ $SPECIES_COLUMNS -eq 0 ]; then
    SPECIES_COLUMNS=$((GCA_COLUMN - 1))
    if [ $SPECIES_COLUMNS -lt 1 ]; then
      echo "Warning: Cannot determine species columns for line: $line"
      continue
    fi
  fi
  
  # Build the species name from the specified number of columns
  SPECIES=""
  for ((i=0; i<SPECIES_COLUMNS; i++)); do
    if [ -n "${FIELDS[$i]}" ]; then
      if [ -n "$SPECIES" ]; then
        SPECIES="$SPECIES ${FIELDS[$i]}"
      else
        SPECIES="${FIELDS[$i]}"
      fi
    fi
  done
  
  # Skip if no species name
  if [ -z "$SPECIES" ]; then
    echo "Warning: Could not extract species name from line: $line"
    continue
  fi
  
  TOTAL=$((TOTAL + 1))
  
  echo "Fetching genome for $SPECIES ($GCA)..."
  
  # Run rapid_fetch.py for this genome
  python MirMachine/bin/rapid_fetch.py -s "$SPECIES" -a "$GCA" -o "$OUTPUT_DIR"
  
  # Check if the command was successful
  if [ $? -eq 0 ]; then
    SUCCESSFUL=$((SUCCESSFUL + 1))
    echo "Successfully fetched genome for $SPECIES ($GCA)"
  else
    echo "Error fetching genome for $SPECIES ($GCA)"
    FAILED=$((FAILED + 1))
  fi
  
done < "$GENOME_FILE"

# Print summary
echo ""
echo "Summary:"
echo "  Total genomes processed: $TOTAL"
echo "  Successfully fetched: $SUCCESSFUL"
echo "  Failed: $FAILED"

exit 0