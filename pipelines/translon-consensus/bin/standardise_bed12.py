import argparse
import sys
import pybedtools as pbt
from Bio.Seq import Seq
from Bio.Data import CodonTable

# --- 1. Constants and Setup ---

STANDARD_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]
STOP_CODONS = STANDARD_TABLE.stop_codons
NEAR_COGNATE_STARTS = STANDARD_TABLE.start_codons + ['GTG', 'TTG', 'ATA', 'ATC', 'ATT'] 

# --- 2. ORF Validation Function ---

def validate_orf(sequence: str) -> bool:
    """
    Validates if a sequence represents a valid Open Reading Frame (ORF).
    Checks: length divisible by 3, start codon, no internal stops, stop codon at end.
    """
    seq_len = len(sequence)
    if seq_len < 3 or seq_len % 3 != 0:
        return False
    
    start_codon = sequence[:3].upper()
    if start_codon not in NEAR_COGNATE_STARTS:
        return False
        
    stop_codon = sequence[-3:].upper()
    if stop_codon not in STOP_CODONS:
        return False
        
    # Check for internal stop codons (must not be present in the coding region)
    internal_seq = sequence[:-3]
    for i in range(0, len(internal_seq), 3):
        codon = internal_seq[i:i+3].upper()
        if codon in STOP_CODONS:
            return False 
            
    return True

# --- 3. Coordinate Conversion Functions ---

def convert_to_0based_slicing(feature: pbt.Interval, is_1based_input: bool):
    """
    Converts a feature's coordinates for sequence extraction if it is 1-based input.
    Returns a new BedTool Interval object ready for 0-based slicing.
    """
    if not is_1based_input:
        # If input is already 0-based, return the original feature
        return feature
    
    # If input is 1-based, we must convert it to 0-based coordinates for slicing
    # 1-based Start is chromStart_input. 0-based Start is chromStart_input - 1
    # 1-based End is chromEnd_input. 0-based End (half-open) is chromEnd_input (same)
    
    new_start = str(feature.start - 1)
    
    # Create a new list of fields with the corrected start
    new_fields = list(feature.fields)
    new_fields[1] = new_start  # Update chromStart (column 2)
    
    # Create a new Interval object from the modified fields for slicing
    return pbt.create_interval_from_list(new_fields)

def convert_to_1based_closed(bed12_line: list) -> list:
    """Converts a standard 0-based, half-open BED12 line to 1-based, fully-closed."""
    try:
        start_0_based = int(bed12_line[1])
        new_start = str(start_0_based + 1)
        new_line = bed12_line[:]
        new_line[1] = new_start
        return new_line
    except ValueError:
        print(f"Error parsing start coordinate in line: {' '.join(bed12_line)}", file=sys.stderr)
        return bed12_line

# --- 4. Main Logic Function with Auto-Detection ---

def process_bed12(bed_path: str, fasta_path: str, output_prefix: str):
    """
    Main function to process the BED12 file, auto-detect coordinate system,
    validate ORFs, and output standardised files.
    """
    
    print(f"Loading BED file: {bed_path}")
    original_bed = pbt.BedTool(bed_path)
    
    # Files to hold the output data
    correct_0_based = []
    correct_1_based = []
    
    print("Auto-detecting input coordinates and validating ORFs...")
    
    for feature in original_bed:
        # ----------------------------------------------------------------------
        # HYPOTHESIS 1: Input is 0-based, half-open (Standard BED)
        # ----------------------------------------------------------------------
        
        # Test 1: Use the feature's coordinates as-is for 0-based slicing
        test_feature_0based = feature 
        
        try:
            # Extract the SPLICED sequence (block=True) for the 0-based interpretation
            seq_0based = test_feature_0based.sequence(fi=fasta_path, s=True, block=True).seq.upper()
            is_valid_0based = validate_orf(seq_0based)
        except Exception:
            is_valid_0based = False

        # ----------------------------------------------------------------------
        # HYPOTHESIS 2: Input is 1-based, fully-closed (Ensembl/GFF style)
        # ----------------------------------------------------------------------
        
        # Test 2: Convert the feature's coordinates to 0-based slicing logic
        test_feature_1based = convert_to_0based_slicing(feature, is_1based_input=True)
        
        try:
            # Extract the SPLICED sequence (block=True) for the 1-based interpretation
            seq_1based = test_feature_1based.sequence(fi=fasta_path, s=True, block=True).seq.upper()
            is_valid_1based = validate_orf(seq_1based)
        except Exception:
            is_valid_1based = False
            
        # ----------------------------------------------------------------------
        # DECISION LOGIC
        # ----------------------------------------------------------------------
        
        # Priority: 0-based > 1-based (since 0-based is the standard BED format)
        
        if is_valid_0based and not is_valid_1based:
            # Input is 0-based and the ORF is correct
            standard_feature = feature
            
        elif is_valid_1based and not is_valid_0based:
            # Input is 1-based, but ORF is correct only after 1-based correction
            standard_feature = test_feature_1based # This feature object is now 0-based!
            print(f"INFO: Converted {feature.name} from 1-based to 0-based.", file=sys.stderr)
            
        elif is_valid_0based and is_valid_1based:
            # Ambiguous: Both interpretations yield a valid ORF (e.g., a very short feature)
            print(f"WARNING: {feature.name} is ambiguous. Defaulting to 0-based.", file=sys.stderr)
            standard_feature = feature
            
        else: # Both failed (not is_valid_0based and not is_valid_1based)
            print(f"Skipping {feature.name}: Failed ORF validation in both systems.", file=sys.stderr)
            continue # Skip to the next feature

        # Feature is now guaranteed to be in the **standard 0-based, half-open** format
        # 1. Output the 0-based standardised feature
        correct_0_based.append(str(standard_feature)) 
        
        # 2. Output the 1-based standardised feature (converted from the 0-based version)
        converted_line = convert_to_1based_closed(standard_feature.fields)
        correct_1_based.append('\t'.join(converted_line))


    # --- 5. Output Writing ---
    
    # Write 0-based file
    output_0_based_path = f"{output_prefix}_0based_standard.bed12"
    with open(output_0_based_path, 'w') as f:
        f.write('\n'.join(correct_0_based) + '\n')
    print(f"\nSuccessfully wrote standardised 0-based BED12 to: {output_0_based_path}")

    # Write 1-based file
    output_1_based_path = f"{output_prefix}_1based_ensembl.bed12"
    with open(output_1_based_path, 'w') as f:
        f.write('\n'.join(correct_1_based) + '\n')
    print(f"Successfully wrote standardised 1-based BED12 to: {output_1_based_path}")

# --- 6. Argument Parsing and Main ---

def main():
    parser = argparse.ArgumentParser(
        description="Auto-detects BED12 coordinate system (0-based vs 1-based) "
                    "by validating the Open Reading Frame (ORF) against a FASTA genome, "
                    "and outputs standardised files."
    )
    # ... (argparse setup remains the same)
    parser.add_argument(
        '-i', '--bed12', 
        required=True, 
        help="Path to the input BED12 file."
    )
    parser.add_argument(
        '-f', '--fasta', 
        required=True, 
        help="Path to the reference genome FASTA file (indexed with samtools faidx)."
    )
    parser.add_argument(
        '-o', '--output_prefix', 
        required=True, 
        help="Prefix for output files."
    )

    args = parser.parse_args()

    # Pybedtools requires a FASTA index (assuming samtools is available)
    try:
        pbt.fasta_idx(args.fasta)
    except Exception as e:
        print(f"Error indexing FASTA file. Error: {e}", file=sys.stderr)
        sys.exit(1)

    process_bed12(args.bed12, args.fasta, args.output_prefix)

if __name__ == '__main__':
    main()