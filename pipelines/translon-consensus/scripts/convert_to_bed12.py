#!/usr/bin/env python3
"""
Auto-detect tool format and convert to BED12.

Detects format based on file structure and calls appropriate converter.
"""

import argparse
import sys
from pathlib import Path
import subprocess


def detect_format(input_file):
    """
    Detect format by examining the file structure.

    Returns: 'ribotie', 'price', 'iribo', 'orfquant', or None
    """
    with open(input_file) as f:
        # Skip comments
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')

            if len(fields) < 9:
                return None

            # Check 2nd column (source) for tool names
            source = fields[1] if len(fields) > 1 else ''

            if source == 'RiboTIE':
                # GTF format: check for gene_id attribute
                if 'gene_id' in fields[-1]:
                    return 'ribotie'

            elif source == 'iRibo':
                # GFF3 format: check for ID= attribute
                if 'ID=' in fields[-1]:
                    return 'iribo'

            elif source == 'ORFQuant':
                # GFF3 format: 9th column is just ORF identifier
                if '=' not in fields[-1] and 'ENST' in fields[-1]:
                    return 'orfquant'

            # Check if it looks like BED12 (PRICE)
            elif len(fields) >= 12:
                # Try to parse as numbers
                try:
                    int(fields[1])  # start
                    int(fields[2])  # end
                    int(fields[9])  # blockCount
                    return 'price'
                except ValueError:
                    pass

            # Only check first data line
            break

    return None


def main():
    parser = argparse.ArgumentParser(
        description='Auto-detect format and convert to BED12',
        epilog="""
Automatically detects tool format (RiboTIE, PRICE, iRibo, ORFQuant)
and converts to proper BED12 format.

Example:
  convert_to_bed12.py input.bed output.bed12
        """
    )
    parser.add_argument('input', help='Input file (any supported format)')
    parser.add_argument('output', help='Output BED12 file')
    parser.add_argument('--tool', choices=['ribotie', 'price', 'iribo', 'orfquant'],
                        help='Force specific tool converter (auto-detect if not specified)')

    args = parser.parse_args()

    # Detect format if not specified
    if args.tool:
        tool = args.tool
        print(f"Using forced tool: {tool}", file=sys.stderr)
    else:
        tool = detect_format(args.input)
        if not tool:
            print("Error: Could not detect file format", file=sys.stderr)
            print("Please specify --tool manually", file=sys.stderr)
            sys.exit(1)
        print(f"Detected format: {tool}", file=sys.stderr)

    # Get script directory
    script_dir = Path(__file__).parent

    # Map tool to converter script
    converters = {
        'ribotie': script_dir / 'ribotie_to_bed12.py',
        'price': script_dir / 'price_to_bed12.py',
        'iribo': script_dir / 'iribo_to_bed12.py',
        'orfquant': script_dir / 'orfquant_to_bed12.py'
    }

    converter_script = converters[tool]

    if not converter_script.exists():
        print(f"Error: Converter script not found: {converter_script}", file=sys.stderr)
        sys.exit(1)

    # Call appropriate converter
    cmd = [str(converter_script), args.input, args.output]
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)

    result = subprocess.run(cmd, capture_output=False)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
