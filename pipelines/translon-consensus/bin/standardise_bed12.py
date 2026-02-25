#!/usr/bin/env python3
"""
Convert various ORF prediction formats to standardized 0-based BED12
and validate ORF sequences.

Handles:
- GFF3/GTF (ORFQuant, iRibo, RiboTIE): 1-based inclusive
- BED12 (PRICE): 0-based half-open [start, end)
"""

import argparse
import sys
from collections import defaultdict, OrderedDict

import pybedtools as pbt


STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG"}  # Canonical start
ALT_START_CODONS = {"CTG", "GTG", "TTG", "ACG", "ATA", "ATT", "ATC"}  # Alternative starts


def parse_args():
    ap = argparse.ArgumentParser(
        description="Convert ORF formats to BED12 and validate sequences"
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Input file (GFF3/GTF/BED)")
    ap.add_argument("-f", "--fasta", required=True,
                    help="Genome FASTA file")
    ap.add_argument("-o", "--output_prefix", required=True,
                    help="Output prefix for BED12 files")
    ap.add_argument("--format", 
                    choices=["gff3", "gtf", "bed12", "auto"],
                    default="auto",
                    help="Input format (default: auto-detect)")
    ap.add_argument("--source", default=None,
                    help="Filter by source/tool name (default: auto-detect from file)")
    ap.add_argument("--min-length", type=int, default=3,
                    help="Minimum ORF length in nucleotides (default: 3)")
    ap.add_argument("--allow-no-stop", action="store_true",
                    help="Allow ORFs without stop codons")
    ap.add_argument("--require-start-codon", action="store_true",
                    help="Require canonical ATG start codon")
    ap.add_argument("--allow-alt-start", action="store_true",
                    help="Allow alternative start codons (CTG, GTG, etc)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def detect_format_and_source(filepath):
    """
    Auto-detect format and source tool from file content.
    Returns: (format, source) tuple
    """
    detected_format = None
    detected_source = None
    
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            
            # Check if BED12 (12 fields, numeric positions)
            if len(fields) == 12:
                try:
                    int(fields[1])
                    int(fields[2])
                    detected_format = "bed12"
                    # Try to infer source from name field patterns
                    name = fields[3]
                    if 'uORF' in name or 'dORF' in name or 'ncRNA' in name:
                        detected_source = "PRICE"
                    break
                except ValueError:
                    pass
            
            # Check if GFF3/GTF (9 fields, feature type in col 3)
            if len(fields) >= 9:
                if fields[2] in ['CDS', 'exon', 'transcript']:
                    source_col = fields[1]
                    
                    # Distinguish GTF from GFF3 by attributes
                    if 'gene_id' in fields[8] and '"' in fields[8]:
                        detected_format = "gtf"
                        # Common GTF tools
                        if source_col == 'RiboTIE' or 'RiboTIE' in fields[8]:
                            detected_source = "RiboTIE"
                        elif source_col == 'Ribo-TISH' or 'TISH' in fields[8]:
                            detected_source = "Ribo-TISH"
                    else:
                        detected_format = "gff3"
                        # Common GFF3 tools
                        if source_col == 'ORFQuant' or 'ORFQuant' in fields[8]:
                            detected_source = "ORFQuant"
                        elif source_col == 'iRibo' or 'iRibo' in fields[8]:
                            detected_source = "iRibo"
                        elif source_col == 'RibORF' or 'RibORF' in fields[8]:
                            detected_source = "RibORF"
                    
                    # Use source column as fallback
                    if not detected_source and source_col != '.':
                        detected_source = source_col
                    
                    break
            
            # Only check first non-comment line
            break
    
    return detected_format, detected_source


def parse_gff3_attributes(attr_string):
    """Parse GFF3 attributes (ID=value;Name=value)"""
    attrs = {}
    
    # Handle case where entire string is just an ID without key=value format
    # e.g., ORFQuant: "ENST00000620552.4_51_6254"
    if '=' not in attr_string:
        # Treat the whole string as an ID
        attrs['ID'] = attr_string.strip()
        return attrs
    
    for item in attr_string.strip().split(';'):
        if not item.strip():
            continue
        if '=' not in item:
            continue
        key, value = item.split('=', 1)
        attrs[key.strip()] = value.strip()
    return attrs


def parse_gtf_attributes(attr_string):
    """Parse GTF attributes (key "value"; key "value";)"""
    attrs = {}
    # Remove trailing semicolon and split
    attr_string = attr_string.strip().rstrip(';')
    
    for item in attr_string.split(';'):
        item = item.strip()
        if not item or ' ' not in item:
            continue
        parts = item.split(None, 1)
        if len(parts) == 2:
            key = parts[0].strip()
            # Strip quotes and handle malformed double quotes
            value = parts[1].strip()
            # Remove all quotes (single and double)
            while value and value[0] in ['"', "'"]:
                value = value[1:]
            while value and value[-1] in ['"', "'"]:
                value = value[:-1]
            attrs[key] = value
    return attrs


def parse_gff_gtf(filepath, file_format, source_filter=None, verbose=False):
    """
    Parse GFF3 or GTF and group exons by ORF/transcript ID.
    Returns: dict of {orf_id: {'chrom', 'strand', 'score', 'exons': [(start, end), ...]}}
    
    Note: GFF3/GTF coordinates are 1-based, inclusive [start, end]
    """
    orfs = defaultdict(lambda: {
        'chrom': None,
        'strand': None,
        'score': '0',
        'source': None,
        'exons': []
    })
    
    line_num = 0
    skipped_lines = 0
    
    with open(filepath) as f:
        for line in f:
            line_num += 1
            
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                if verbose:
                    print(f"[WARN] Line {line_num}: Expected 9 fields, got {len(fields)}", 
                          file=sys.stderr)
                skipped_lines += 1
                continue
            
            # Strip whitespace from all fields (handles RiboTIE's space-padded coordinates)
            fields = [f.strip() for f in fields]
            
            chrom = fields[0]
            source = fields[1]
            feature = fields[2]
            
            try:
                start = int(fields[3])  # 1-based
                end = int(fields[4])    # 1-based, inclusive
            except ValueError as e:
                if verbose:
                    print(f"[WARN] Line {line_num}: Invalid coordinates: {e}", file=sys.stderr)
                skipped_lines += 1
                continue
            
            score = fields[5]
            strand = fields[6]
            attr_string = fields[8]
            
            # Filter by source if requested
            if source_filter and source != source_filter:
                continue
            
            # Only process CDS features
            if feature != 'CDS':
                continue
            
            # Parse attributes to get ORF ID
            if file_format == "gff3":
                attrs = parse_gff3_attributes(attr_string)
                orf_id = attrs.get('ID', attrs.get('Parent', f'unknown_{chrom}_{start}'))
            else:  # GTF
                attrs = parse_gtf_attributes(attr_string)
                orf_id = attrs.get('ORF_id', attrs.get('transcript_id', f'unknown_{chrom}_{start}'))
            
            # Convert to 0-based half-open [start, end)
            start_0based = start - 1
            end_0based = end  # End is already exclusive in half-open notation
            
            # Store exon
            orf_data = orfs[orf_id]
            orf_data['chrom'] = chrom
            orf_data['strand'] = strand
            orf_data['source'] = source
            if score != '.':
                orf_data['score'] = score
            orf_data['exons'].append((start_0based, end_0based))
            
            if verbose:
                print(f"[DEBUG] {orf_id}: exon {start}-{end} (1-based) -> "
                      f"{start_0based}-{end_0based} (0-based)", file=sys.stderr)
    
    if skipped_lines > 0:
        print(f"[WARN] Skipped {skipped_lines} malformed lines", file=sys.stderr)
    
    return orfs


def parse_bed12(filepath, verbose=False):
    """
    Parse BED12 file. Assumes already 0-based half-open.
    Returns same structure as parse_gff_gtf for consistency.
    
    Also fixes common BED12 issues:
    - thickStart == thickEnd == chromStart (no CDS) -> set thickEnd = chromEnd
    """
    orfs = {}
    line_num = 0
    fixed_thick = 0
    
    with open(filepath) as f:
        for line in f:
            line_num += 1
            
            if line.startswith('#') or line.startswith('track') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 12:
                if verbose:
                    print(f"[WARN] Line {line_num}: Expected 12 fields, got {len(fields)}", 
                          file=sys.stderr)
                continue
            
            chrom = fields[0]
            start = int(fields[1])  # Already 0-based
            end = int(fields[2])    # Already 0-based
            name = fields[3]
            score = fields[4]
            strand = fields[5]
            thick_start = int(fields[6])
            thick_end = int(fields[7])
            block_count = int(fields[9])
            block_sizes = [int(x) for x in fields[10].rstrip(',').split(',') if x]
            block_starts = [int(x) for x in fields[11].rstrip(',').split(',') if x]
            
            # Fix common issue: thickStart == thickEnd == chromStart (indicates no CDS)
            # For ORFs, we want the entire feature to be thick (coding)
            if thick_start == thick_end == start:
                thick_end = end
                fixed_thick += 1
                if verbose:
                    print(f"[FIX] {name}: Fixed thickEnd {start} -> {end}", file=sys.stderr)
            
            # Reconstruct exons
            exons = []
            for i in range(block_count):
                exon_start = start + block_starts[i]
                exon_end = exon_start + block_sizes[i]
                exons.append((exon_start, exon_end))
            
            orfs[name] = {
                'chrom': chrom,
                'strand': strand,
                'score': score,
                'source': 'bed12',
                'exons': exons,
                'thick_start': thick_start,
                'thick_end': thick_end
            }
            
            if verbose:
                print(f"[DEBUG] {name}: {len(exons)} exons", file=sys.stderr)
    
    if fixed_thick > 0:
        print(f"[INFO] Fixed thickEnd for {fixed_thick} ORFs (was thickStart==thickEnd==chromStart)", 
              file=sys.stderr)
    
    return orfs


def validate_bed12_line(line, line_num=0):
    """
    Validate that a BED12 line has proper format.
    Returns: (is_valid, error_message)
    """
    fields = line.strip().split('\t')
    
    if len(fields) != 12:
        return False, f"Expected 12 fields, got {len(fields)}"
    
    # Check required numeric fields
    try:
        chrom_start = int(fields[1])
        chrom_end = int(fields[2])
        thick_start = int(fields[6])
        thick_end = int(fields[7])
        block_count = int(fields[9])
        
        if chrom_start < 0 or chrom_end < 0:
            return False, "Negative coordinates"
        
        if chrom_start >= chrom_end:
            return False, f"Invalid range: start={chrom_start} >= end={chrom_end}"
        
        if block_count < 1:
            return False, f"Invalid block_count: {block_count}"
        
    except (ValueError, IndexError) as e:
        return False, f"Invalid numeric field: {e}"
    
    # Check block information
    try:
        block_sizes = [int(x) for x in fields[10].rstrip(',').split(',') if x]
        block_starts = [int(x) for x in fields[11].rstrip(',').split(',') if x]
        
        if len(block_sizes) != block_count:
            return False, f"blockSizes count {len(block_sizes)} != blockCount {block_count}"
        
        if len(block_starts) != block_count:
            return False, f"blockStarts count {len(block_starts)} != blockCount {block_count}"
        
        if any(s <= 0 for s in block_sizes):
            return False, "Block sizes must be positive"
        
        if block_starts[0] != 0:
            return False, f"First blockStart must be 0, got {block_starts[0]}"
        
        # Check that blocks don't overlap and are within bounds
        for i in range(block_count):
            block_end = block_starts[i] + block_sizes[i]
            if block_end > (chrom_end - chrom_start):
                return False, f"Block {i} extends beyond feature end"
            
            if i < block_count - 1:
                if block_starts[i+1] < block_end:
                    return False, f"Blocks {i} and {i+1} overlap"
        
    except (ValueError, IndexError) as e:
        return False, f"Invalid block information: {e}"
    
    # Check strand
    if fields[5] not in ['+', '-', '.']:
        return False, f"Invalid strand: {fields[5]}"
    
    return True, "OK"


def orfs_to_bed12(orfs, chrom_mapping=None, verbose=False):
    """
    Convert ORF structure to BED12 format.
    Returns: dict of {orf_id: bed12_line}
    
    Args:
        orfs: dict of ORF data
        chrom_mapping: dict mapping input chrom names to FASTA chrom names
        verbose: print debug info
    """
    bed12_dict = OrderedDict()
    unmapped_chroms = set()
    invalid_orfs = []
    
    for orf_id, data in orfs.items():
        exons = sorted(data['exons'], key=lambda x: x[0])
        
        if not exons:
            invalid_orfs.append((orf_id, "no_exons"))
            continue
        
        chrom = data['chrom']
        
        # Apply chromosome name mapping if provided
        if chrom_mapping:
            if chrom in chrom_mapping:
                chrom = chrom_mapping[chrom]
            else:
                unmapped_chroms.add(data['chrom'])
                if verbose:
                    print(f"[WARN] {orf_id}: chromosome '{data['chrom']}' not found in FASTA", 
                          file=sys.stderr)
                invalid_orfs.append((orf_id, f"unmapped_chrom:{data['chrom']}"))
                continue
        
        strand = data['strand']
        score = data['score']
        
        # Calculate transcript boundaries
        tx_start = exons[0][0]
        tx_end = exons[-1][1]
        
        # For ORFs, thick_start = tx_start, thick_end = tx_end (full CDS)
        thick_start = tx_start
        thick_end = tx_end
        
        # RGB color (default black)
        rgb = "0,0,0"
        
        # Block information
        block_count = len(exons)
        block_sizes = [str(end - start) for start, end in exons]
        block_starts = [str(start - tx_start) for start, end in exons]
        
        # Build BED12 line
        bed12_line = '\t'.join([
            chrom,
            str(tx_start),
            str(tx_end),
            orf_id,
            score,
            strand,
            str(thick_start),
            str(thick_end),
            rgb,
            str(block_count),
            ','.join(block_sizes) + ',',
            ','.join(block_starts) + ','
        ])
        
        # Validate the generated BED12 line
        is_valid, error = validate_bed12_line(bed12_line)
        if not is_valid:
            invalid_orfs.append((orf_id, f"invalid_bed12:{error}"))
            if verbose:
                print(f"[ERROR] {orf_id}: Generated invalid BED12 - {error}", file=sys.stderr)
            continue
        
        bed12_dict[orf_id] = bed12_line
        
        if verbose:
            total_len = sum(int(s) for s in block_sizes)
            mapped_note = f" (mapped from {data['chrom']})" if chrom_mapping and data['chrom'] != chrom else ""
            print(f"[INFO] {orf_id}: {chrom}:{tx_start}-{tx_end}({strand}) "
                  f"{block_count} exons, {total_len}bp{mapped_note}", file=sys.stderr)
    
    # Report issues
    if invalid_orfs:
        print(f"[WARN] Skipped {len(invalid_orfs)} invalid ORFs", file=sys.stderr)
        if verbose:
            for orf_id, reason in invalid_orfs[:10]:
                print(f"  {orf_id}: {reason}", file=sys.stderr)
            if len(invalid_orfs) > 10:
                print(f"  ... and {len(invalid_orfs)-10} more", file=sys.stderr)
    
    if unmapped_chroms and not verbose:
        print(f"[WARN] {len(unmapped_chroms)} chromosomes not found in FASTA: "
              f"{', '.join(sorted(list(unmapped_chroms)[:5]))}"
              f"{' ...' if len(unmapped_chroms) > 5 else ''}", file=sys.stderr)
    
    return bed12_dict


def load_fasta_sequences(fasta_path):
    """Load FASTA into {header: sequence}"""
    seqs = OrderedDict()
    name = None
    buf = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = line[1:]
                buf = []
            else:
                buf.append(line)

        if name:
            seqs[name] = "".join(buf).upper()

    return seqs


def get_fasta_chromosomes(fasta_path):
    """
    Extract chromosome names from FASTA file.
    Returns: set of chromosome names
    """
    chroms = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                # Take first word after '>' as chrom name
                chrom = line[1:].split()[0]
                chroms.add(chrom)
    return chroms


def build_chrom_mapping(input_chroms, fasta_chroms):
    """
    Build mapping between input chromosome names and FASTA chromosome names.
    Handles common naming variations:
    - chr1 <-> 1
    - chrM <-> MT
    - chr_scaffold <-> scaffold
    
    Returns: dict mapping input_chrom -> fasta_chrom
    """
    mapping = {}
    
    # First, try exact matches
    for chrom in input_chroms:
        if chrom in fasta_chroms:
            mapping[chrom] = chrom
    
    # For unmatched chroms, try common conversions
    unmatched = input_chroms - set(mapping.keys())
    available = fasta_chroms - set(mapping.values())
    
    for chrom in unmatched:
        # Try adding 'chr' prefix
        if chrom not in fasta_chroms and f"chr{chrom}" in available:
            mapping[chrom] = f"chr{chrom}"
            available.remove(f"chr{chrom}")
            continue
        
        # Try removing 'chr' prefix
        if chrom.startswith('chr'):
            no_chr = chrom[3:]
            if no_chr in available:
                mapping[chrom] = no_chr
                available.remove(no_chr)
                continue
        
        # Special case: chrM <-> MT, chrMT <-> MT
        if chrom in ['chrM', 'M']:
            if 'MT' in available:
                mapping[chrom] = 'MT'
                available.remove('MT')
                continue
            elif 'chrMT' in available:
                mapping[chrom] = 'chrMT'
                available.remove('chrMT')
                continue
        elif chrom in ['chrMT', 'MT']:
            if 'M' in available:
                mapping[chrom] = 'M'
                available.remove('M')
                continue
            elif 'chrM' in available:
                mapping[chrom] = 'chrM'
                available.remove('chrM')
                continue
    
    return mapping


def extract_spliced_sequences(bed12_lines, fasta_path):
    """
    Extract sequences for BED12 lines using bedtools.
    Returns: dict[orf_id -> sequence]
    """
    import tempfile
    import os
    
    # Write BED12 to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
        tmp_path = tmp.name
        for line in bed12_lines.values():
            tmp.write(line + '\n')
    
    try:
        bt = pbt.BedTool(tmp_path)
        fa = bt.sequence(fi=fasta_path, s=True, split=True, name=True)
        seqs = load_fasta_sequences(fa.seqfn)
        
        # Map back to ORF IDs (bedtools adds ::chrom:start-end(strand))
        orf_seqs = {}
        for header, seq in seqs.items():
            # Extract orf_id from "orf_id::chrom:start-end(strand)"
            orf_id = header.split('::')[0] if '::' in header else header
            orf_seqs[orf_id] = seq
        
        return orf_seqs
    finally:
        os.unlink(tmp_path)


def get_downstream_sequence(chrom, end, strand, fasta_path, length=3):
    """
    Get downstream sequence after the ORF end.
    For + strand: sequence from end to end+length
    For - strand: sequence from end-length to end (reverse complement)
    
    Returns: sequence string or None if error
    """
    import tempfile
    import os
    
    # Create a BED entry for the downstream region
    if strand == '+':
        # For + strand, get bases after the end
        bed_line = f"{chrom}\t{end}\t{end+length}\tdownstream\t0\t{strand}\n"
    else:
        # For - strand, get bases before the start (which is stored as 'end' in 0-based)
        # Actually, 'end' is the chromEnd, so for - strand we want bases before chromStart
        # Wait, we need chromStart for this. Let me redesign...
        return None
    
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
            tmp_path = tmp.name
            tmp.write(bed_line)
        
        bt = pbt.BedTool(tmp_path)
        fa = bt.sequence(fi=fasta_path, s=True)
        seqs = load_fasta_sequences(fa.seqfn)
        
        os.unlink(tmp_path)
        
        if seqs:
            return list(seqs.values())[0]
        return None
    except:
        return None


def extend_orf_for_stop_codon(orf_id, orf_data, fasta_path, verbose=False):
    """
    Check multiple possibilities for stop codon location and extend if found:
    1. Already included in coordinates (last 3bp)
    2. In next 3bp after end (tool excluded stop)
    3. Overlapping by 1-2bp (coordinate system issue)
    
    This handles uncertainties in whether coordinates are 0-based or 1-based,
    and whether the stop codon is included or excluded.
    
    Returns: (extended, reason) where extended is bool
    """
    exons = orf_data['exons']
    strand = orf_data['strand']
    chrom = orf_data['chrom']
    
    if not exons:
        return False, "no_exons"
    
    # Sort exons by position
    sorted_exons = sorted(exons, key=lambda x: x[0])
    
    import tempfile
    import os
    
    # Define regions to check based on strand
    if strand == '+':
        # For + strand, check around the end of the last exon
        last_exon_idx = len(sorted_exons) - 1
        last_exon_start, last_exon_end = sorted_exons[last_exon_idx]
        
        # Check multiple windows:
        # 1. Last 3bp of current ORF (stop already included?)
        # 2. Next 3bp after ORF (stop excluded?)
        # 3. Next 1-5bp (off-by-one errors)
        check_regions = [
            (last_exon_end - 3, last_exon_end, 0, "already_included"),
            (last_exon_end, last_exon_end + 3, 3, "extend_3bp"),
            (last_exon_end + 1, last_exon_end + 4, 4, "extend_4bp_off1"),
            (last_exon_end + 2, last_exon_end + 5, 5, "extend_5bp_off2"),
            (last_exon_end - 1, last_exon_end + 2, 2, "extend_2bp_overlap1"),
        ]
    else:  # - strand
        # For - strand, check around the start of the first exon
        last_exon_idx = 0
        last_exon_start, last_exon_end = sorted_exons[last_exon_idx]
        
        # Check multiple windows (going backwards on - strand)
        check_regions = [
            (last_exon_start, last_exon_start + 3, 0, "already_included"),
            (max(0, last_exon_start - 3), last_exon_start, 3, "extend_3bp"),
            (max(0, last_exon_start - 4), last_exon_start - 1, 4, "extend_4bp_off1"),
            (max(0, last_exon_start - 5), last_exon_start - 2, 5, "extend_5bp_off2"),
            (max(0, last_exon_start - 2), last_exon_start + 1, 2, "extend_2bp_overlap1"),
        ]
    
    # Try each region
    for check_start, check_end, extend_by, label in check_regions:
        try:
            bed_line = f"{chrom}\t{check_start}\t{check_end}\ttemp\t0\t{strand}"
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
                tmp_path = tmp.name
                tmp.write(bed_line + '\n')
            
            bt = pbt.BedTool(tmp_path)
            fa = bt.sequence(fi=fasta_path, s=True)
            seqs = load_fasta_sequences(fa.seqfn)
            os.unlink(tmp_path)
            
            if not seqs:
                continue
            
            seq = list(seqs.values())[0].upper()
            
            # Check if it's exactly 3bp and is a stop codon
            if len(seq) == 3 and seq in STOP_CODONS:
                # Found a stop codon!
                if extend_by > 0:
                    # Need to extend the ORF
                    if strand == '+':
                        sorted_exons[last_exon_idx] = (last_exon_start, last_exon_end + extend_by)
                    else:
                        sorted_exons[last_exon_idx] = (last_exon_start - extend_by, last_exon_end)
                    
                    orf_data['exons'] = sorted_exons
                    
                    if verbose:
                        print(f"[EXTEND] {orf_id}: {label} - found {seq} at "
                              f"{chrom}:{check_start}-{check_end}({strand})", file=sys.stderr)
                    
                    return True, f"{label}_{seq}"
                else:
                    # Stop already included
                    if verbose:
                        print(f"[OK] {orf_id}: Stop codon {seq} already included", file=sys.stderr)
                    return False, f"stop_included_{seq}"
        
        except Exception as e:
            if verbose:
                print(f"[DEBUG] {orf_id}: Error checking {label}: {e}", file=sys.stderr)
            continue
    
    # No stop codon found in any position
    return False, "no_stop_found"


def has_internal_stops(seq):
    """Check for internal stop codons (excluding the last codon)"""
    num_codons = len(seq) // 3
    
    for i in range(0, (num_codons - 1) * 3, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            return True, i
    
    return False, None


def validate_orf(seq, min_length=3, allow_no_stop=False, require_start=False, allow_alt_start=False):
    """
    Strict ORF validation:
      - length >= min_length
      - length multiple of 3
      - no internal stop codons
      - terminal stop codon (unless allow_no_stop)
      - start codon (if require_start)
    
    Returns: (is_valid, reason, details)
    """
    seq_len = len(seq)
    
    if seq_len < min_length:
        return False, "too_short", f"len={seq_len}"
    
    if seq_len % 3 != 0:
        return False, "not_multiple_of_3", f"len={seq_len},mod={seq_len%3}"
    
    # Check start codon if required
    if require_start:
        first_codon = seq[:3]
        has_canonical_start = first_codon in START_CODONS
        has_alt_start = first_codon in ALT_START_CODONS
        
        if not has_canonical_start:
            if allow_alt_start and has_alt_start:
                # Alternative start is OK
                pass
            else:
                return False, "no_valid_start", f"first_codon={first_codon}"
    
    # Check for internal stops
    has_internal, pos = has_internal_stops(seq)
    if has_internal:
        codon = seq[pos:pos+3]
        return False, "internal_stop", f"pos={pos},codon={codon}"
    
    # Check terminal stop codon
    last_codon = seq[-3:]
    has_stop = last_codon in STOP_CODONS
    
    if not has_stop and not allow_no_stop:
        return False, "no_terminal_stop", f"last_codon={last_codon}"
    
    # Determine validation type
    first_codon = seq[:3]
    if has_stop:
        if first_codon in START_CODONS:
            return True, "valid_canonical", f"len={seq_len},start=ATG,stop={last_codon}"
        elif first_codon in ALT_START_CODONS:
            return True, "valid_alt_start", f"len={seq_len},start={first_codon},stop={last_codon}"
        else:
            return True, "valid_no_start", f"len={seq_len},first={first_codon},stop={last_codon}"
    else:
        return True, "valid_no_stop", f"len={seq_len}"


def fix_frameshift_issues(orf_id, orf_data, orf_seqs, require_start=False, allow_alt_start=False, verbose=False):
    """
    If an ORF fails validation due to not being multiple of 3,
    try trimming 1-2bp from start or end to find a valid frame.
    
    Also checks if the issue is at exon boundaries (splice junctions).
    
    Modifies orf_data in place if a fix is found.
    Returns: (fixed, adjustment_made)
    """
    seq = orf_seqs.get(orf_id)
    if not seq:
        return False, "no_seq"
    
    # First check if it's already valid
    is_valid, reason, details = validate_orf(seq, 3, False, require_start, allow_alt_start)
    if is_valid:
        return False, "already_valid"
    
    # Only try to fix frameshift issues
    if "not_multiple_of_3" not in reason:
        return False, reason
    
    seq_len = len(seq)
    remainder = seq_len % 3
    exons = sorted(orf_data['exons'], key=lambda x: x[0])
    strand = orf_data['strand']
    num_exons = len(exons)
    
    # Check if this is a multi-exon ORF (potential splice junction issue)
    at_junction = num_exons > 1
    
    # Try different adjustments
    adjustments = []
    
    # If length % 3 == 1, try removing 1bp from start or end
    if remainder == 1:
        adjustments = [
            (1, 0, f"trim_1_from_start{'_junction' if at_junction else ''}"),
            (0, 1, f"trim_1_from_end{'_junction' if at_junction else ''}"),
        ]
    # If length % 3 == 2, try removing 2bp from start/end, or 1bp from each
    elif remainder == 2:
        adjustments = [
            (2, 0, f"trim_2_from_start{'_junction' if at_junction else ''}"),
            (0, 2, f"trim_2_from_end{'_junction' if at_junction else ''}"),
            (1, 1, f"trim_1_from_each{'_junction' if at_junction else ''}"),
        ]
    
    for trim_start, trim_end, label in adjustments:
        test_seq = seq[trim_start:] if trim_end == 0 else seq[trim_start:-trim_end]
        
        # Validate the adjusted sequence
        is_valid, val_reason, val_details = validate_orf(test_seq, 3, False, require_start, allow_alt_start)
        
        if is_valid:
            # Found a valid frame! Adjust the ORF coordinates
            if strand == '+':
                # Trim from genomic start (first exon)
                if trim_start > 0:
                    first_start, first_end = exons[0]
                    exons[0] = (first_start + trim_start, first_end)
                
                # Trim from genomic end (last exon)
                if trim_end > 0:
                    last_start, last_end = exons[-1]
                    exons[-1] = (last_start, last_end - trim_end)
            
            else:  # - strand (coordinates are reversed)
                # For - strand: trim_start affects the last exon (3' end in transcript)
                if trim_start > 0:
                    last_start, last_end = exons[-1]
                    exons[-1] = (last_start, last_end - trim_start)
                
                # trim_end affects the first exon (5' end in transcript)
                if trim_end > 0:
                    first_start, first_end = exons[0]
                    exons[0] = (first_start + trim_end, first_end)
            
            orf_data['exons'] = exons
            
            if verbose:
                junction_note = f" (multi-exon: {num_exons} exons)" if at_junction else ""
                print(f"[FIX] {orf_id}: {label}{junction_note} - {seq_len}bp -> {len(test_seq)}bp, "
                      f"now {val_reason}", file=sys.stderr)
            
            return True, label
    
    return False, f"no_fix_found_mod{remainder}{'_junction' if at_junction else ''}"


def main():
    args = parse_args()

    # Auto-detect format and source
    file_format = args.format
    source_filter = args.source
    
    if file_format == "auto" or source_filter is None:
        detected_format, detected_source = detect_format_and_source(args.input)
        
        if file_format == "auto":
            if not detected_format:
                print("[ERROR] Could not auto-detect format", file=sys.stderr)
                sys.exit(1)
            file_format = detected_format
            print(f"[INFO] Detected format: {file_format}", file=sys.stderr)
        
        if source_filter is None and detected_source:
            source_filter = detected_source
            print(f"[INFO] Detected source: {source_filter}", file=sys.stderr)
    
    # Step 1: Parse input and convert to BED12
    print(f"[INFO] Reading {args.input} as {file_format}", file=sys.stderr)
    if source_filter:
        print(f"[INFO] Filtering for source: {source_filter}", file=sys.stderr)
    
    if file_format in ["gff3", "gtf"]:
        orfs = parse_gff_gtf(args.input, file_format, source_filter, args.verbose)
    elif file_format == "bed12":
        orfs = parse_bed12(args.input, args.verbose)
    else:
        print(f"[ERROR] Unsupported format: {file_format}", file=sys.stderr)
        sys.exit(1)
    
    print(f"[INFO] Parsed {len(orfs)} ORFs", file=sys.stderr)
    
    # Get chromosome names from input
    input_chroms = {data['chrom'] for data in orfs.values() if data['chrom']}
    print(f"[INFO] Input chromosomes: {len(input_chroms)} unique", file=sys.stderr)
    
    # Get chromosome names from FASTA
    print(f"[INFO] Reading chromosome names from {args.fasta}", file=sys.stderr)
    fasta_chroms = get_fasta_chromosomes(args.fasta)
    print(f"[INFO] FASTA chromosomes: {len(fasta_chroms)} unique", file=sys.stderr)
    
    # Build chromosome mapping
    chrom_mapping = build_chrom_mapping(input_chroms, fasta_chroms)
    
    exact_matches = sum(1 for k, v in chrom_mapping.items() if k == v)
    converted_matches = len(chrom_mapping) - exact_matches
    unmapped = len(input_chroms) - len(chrom_mapping)
    
    print(f"[INFO] Chromosome mapping: {exact_matches} exact, {converted_matches} converted, "
          f"{unmapped} unmapped", file=sys.stderr)
    
    if args.verbose and converted_matches > 0:
        print("[INFO] Chromosome conversions:", file=sys.stderr)
        for in_chr, out_chr in chrom_mapping.items():
            if in_chr != out_chr:
                print(f"  {in_chr} -> {out_chr}", file=sys.stderr)
    
    if unmapped > 0:
        unmapped_list = sorted(input_chroms - set(chrom_mapping.keys()))
        print(f"[WARN] Unmapped chromosomes: {', '.join(unmapped_list[:10])}"
              f"{' ...' if len(unmapped_list) > 10 else ''}", file=sys.stderr)
    
    # Apply chromosome mapping to all ORFs FIRST
    orfs_to_remove = []
    for orf_id, orf_data in orfs.items():
        if chrom_mapping:
            original_chrom = orf_data['chrom']
            if original_chrom in chrom_mapping:
                orf_data['chrom'] = chrom_mapping[original_chrom]
            else:
                # Mark for removal - chromosome not in FASTA
                orfs_to_remove.append(orf_id)
    
    # Remove ORFs on unmapped chromosomes
    for orf_id in orfs_to_remove:
        del orfs[orf_id]
    
    if orfs_to_remove:
        print(f"[INFO] Removed {len(orfs_to_remove)} ORFs on unmapped chromosomes", file=sys.stderr)
    
    # Step 1.5: Check for missing stop codons and extend coordinates if needed
    print(f"[INFO] Checking for stop codons and extending coordinates if needed", file=sys.stderr)
    extended_count = 0
    extension_reasons = {}
    
    for orf_id, orf_data in orfs.items():
        extended, reason = extend_orf_for_stop_codon(orf_id, orf_data, args.fasta, args.verbose)
        if extended:
            extended_count += 1
        extension_reasons[reason] = extension_reasons.get(reason, 0) + 1
    
    print(f"[INFO] Extended {extended_count} ORFs to include stop codons", file=sys.stderr)
    if extension_reasons:
        print("[INFO] Extension results (top 10):", file=sys.stderr)
        for reason, count in sorted(extension_reasons.items(), key=lambda x: -x[1])[:10]:
            print(f"  {reason}: {count}", file=sys.stderr)
    
    # Convert to BED12 format (chromosome names already mapped)
    bed12_dict = orfs_to_bed12(orfs, None, args.verbose)  # Pass None since already mapped
    
    # Step 2: Extract sequences using bedtools
    print(f"[INFO] Extracting spliced sequences from {args.fasta}", file=sys.stderr)
    orf_seqs = extract_spliced_sequences(bed12_dict, args.fasta)
    print(f"[INFO] Extracted {len(orf_seqs)} sequences", file=sys.stderr)
    
    # Step 2.5: Fix frameshift issues by trimming 1-2bp from start/end
    print(f"[INFO] Checking for frameshift issues and adjusting coordinates", file=sys.stderr)
    fixed_count = 0
    fix_reasons = {}
    
    for orf_id in list(bed12_dict.keys()):
        if orf_id not in orfs:
            continue
        
        fixed, reason = fix_frameshift_issues(orf_id, orfs[orf_id], orf_seqs, 
                                              args.require_start_codon, args.allow_alt_start, 
                                              args.verbose)
        if fixed:
            fixed_count += 1
            # Regenerate BED12 line with adjusted coordinates
            bed12_dict = orfs_to_bed12(orfs, None, args.verbose)
        fix_reasons[reason] = fix_reasons.get(reason, 0) + 1
    
    print(f"[INFO] Fixed {fixed_count} ORFs with frameshift issues", file=sys.stderr)
    if fix_reasons:
        print("[INFO] Frameshift fix results (top 10):", file=sys.stderr)
        for reason, count in sorted(fix_reasons.items(), key=lambda x: -x[1])[:10]:
            print(f"  {reason}: {count}", file=sys.stderr)
    
    # Re-extract sequences after fixes
    if fixed_count > 0:
        print(f"[INFO] Re-extracting sequences after frameshift fixes", file=sys.stderr)
        orf_seqs = extract_spliced_sequences(bed12_dict, args.fasta)
    
    # Step 3: Validate sequences
    print(f"[INFO] Validating ORF sequences (min_length={args.min_length}bp)", 
          file=sys.stderr)
    
    valid_orfs = []
    invalid_orfs = []
    validation_counts = {}
    rejection_counts = {}
    
    for orf_id, bed12_line in bed12_dict.items():
        seq = orf_seqs.get(orf_id)
        
        if seq is None:
            rejection_counts["no_sequence"] = rejection_counts.get("no_sequence", 0) + 1
            invalid_orfs.append((orf_id, bed12_line, "no_sequence", ""))
            if args.verbose:
                print(f"[INVALID] {orf_id}: no sequence found", file=sys.stderr)
            continue
        
        is_valid, reason, details = validate_orf(seq, args.min_length, args.allow_no_stop,
                                                  args.require_start_codon, args.allow_alt_start)
        
        if is_valid:
            validation_counts[reason] = validation_counts.get(reason, 0) + 1
            valid_orfs.append((orf_id, bed12_line))
            
            if args.verbose:
                print(f"[VALID] {orf_id}: {reason} ({details})", file=sys.stderr)
        else:
            rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
            invalid_orfs.append((orf_id, bed12_line, reason, details))
            
            if args.verbose:
                print(f"[INVALID] {orf_id}: {reason} ({details})", file=sys.stderr)

    # Step 3.5: Deduplicate by genomic coordinates
    print(f"[INFO] Deduplicating ORFs by genomic coordinates", file=sys.stderr)

    # Parse BED12 lines to extract coordinates
    # BED12 format: chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
    coord_to_orfs = defaultdict(list)

    for orf_id, bed12_line in valid_orfs:
        fields = bed12_line.strip().split('\t')
        if len(fields) >= 12:
            # Use chr, thickStart, thickEnd, strand, blockCount, blockSizes, blockStarts as key
            # This uniquely identifies the genomic ORF structure
            chrom = fields[0]
            thick_start = fields[6]
            thick_end = fields[7]
            strand = fields[5]
            block_count = fields[9]
            block_sizes = fields[10]
            block_starts = fields[11]

            coord_key = (chrom, thick_start, thick_end, strand, block_count, block_sizes, block_starts)
            coord_to_orfs[coord_key].append((orf_id, bed12_line))

    # Keep one representative per unique genomic coordinate
    dedup_valid_orfs = []
    duplicates_removed = 0

    for coord_key, orf_list in coord_to_orfs.items():
        if len(orf_list) > 1:
            duplicates_removed += len(orf_list) - 1
            if args.verbose:
                orf_ids = [orf_id for orf_id, _ in orf_list]
                print(f"[DEDUP] {len(orf_list)} identical ORFs at {coord_key[0]}:{coord_key[1]}-{coord_key[2]}: {', '.join(orf_ids)}", file=sys.stderr)

        # Keep the first one (arbitrary choice)
        dedup_valid_orfs.append(orf_list[0])

    print(f"[INFO] Removed {duplicates_removed} duplicate ORFs ({len(valid_orfs)} → {len(dedup_valid_orfs)})", file=sys.stderr)

    # Replace valid_orfs with deduplicated version
    valid_orfs = dedup_valid_orfs

    # Step 4: Write output files
    valid_file = f"{args.output_prefix}.valid.bed12"
    invalid_file = f"{args.output_prefix}.invalid.bed12"
    
    with open(valid_file, 'w') as out:
        for orf_id, bed12_line in valid_orfs:
            out.write(bed12_line + '\n')
    
    with open(invalid_file, 'w') as out:
        for orf_id, bed12_line, reason, details in invalid_orfs:
            out.write(bed12_line + '\n')
    
    # Step 5: Validate output file integrity
    print(f"\n[INFO] Validating output files...", file=sys.stderr)
    
    valid_line_count = 0
    valid_errors = []
    with open(valid_file, 'r') as f:
        for i, line in enumerate(f, 1):
            is_valid, error = validate_bed12_line(line, i)
            if not is_valid:
                valid_errors.append((i, error))
            else:
                valid_line_count += 1
    
    if valid_errors:
        print(f"[ERROR] Found {len(valid_errors)} invalid lines in {valid_file}:", file=sys.stderr)
        for line_num, error in valid_errors[:5]:
            print(f"  Line {line_num}: {error}", file=sys.stderr)
        if len(valid_errors) > 5:
            print(f"  ... and {len(valid_errors)-5} more errors", file=sys.stderr)
    else:
        print(f"[OK] {valid_file}: All {valid_line_count} lines are valid BED12", file=sys.stderr)
    
    # Step 6: Print summary
    total = len(bed12_dict)
    kept = len(valid_orfs)
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"VALIDATION SUMMARY", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Total ORFs: {total}", file=sys.stderr)
    print(f"Valid ORFs: {kept} ({100*kept/total:.1f}%)", file=sys.stderr)
    print(f"Invalid ORFs: {total-kept} ({100*(total-kept)/total:.1f}%)", file=sys.stderr)
    
    if validation_counts:
        print(f"\nValidation types:", file=sys.stderr)
        for vtype, count in sorted(validation_counts.items(), key=lambda x: -x[1]):
            pct = 100 * count / kept if kept > 0 else 0
            print(f"  {vtype}: {count} ({pct:.1f}%)", file=sys.stderr)
    
    if rejection_counts:
        print(f"\nRejection reasons:", file=sys.stderr)
        for reason, count in sorted(rejection_counts.items(), key=lambda x: -x[1]):
            pct = 100 * count / total
            print(f"  {reason}: {count} ({pct:.1f}%)", file=sys.stderr)
    
    print(f"\nOutput files:", file=sys.stderr)
    print(f"  Valid: {valid_file}", file=sys.stderr)
    print(f"  Invalid: {invalid_file}", file=sys.stderr)
    if kept < (total * 0.1):
        raise Exception(f"insufficient kept ({kept}/{total}) - check logs for issues")


if __name__ == "__main__":
    main()