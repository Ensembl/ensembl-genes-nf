#!/usr/bin/env python3
"""
Process BAM files using length-specific offsets from a tab-delimited file
Optimized version using numpy arrays and buffered writes
"""

import argparse
import pysam 
import numpy as np
from pathlib import Path
import logging
from typing import Dict, Tuple
import itertools
from itertools import groupby
from operator import itemgetter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Default values for configuration
DEFAULT_CHUNK_SIZE = 100000  # Size of chunks to process at once
DEFAULT_WRITE_BUFFER = 10000  # Number of positions to buffer before writing

# Global variables that can be modified via command line
chunk_size = DEFAULT_CHUNK_SIZE
write_buffer = DEFAULT_WRITE_BUFFER

def read_offset_file(offset_file: str) -> Dict[int, int]:
    '''Read offset file and return dictionary of read length to offset'''
    offsets = {}
    with open(offset_file) as f:
        next(f)  # Skip header
        for line in f:
            read_length, offset = line.strip().split('\t')
            offsets[int(read_length)] = int(offset)
    return offsets

class CoverageAccumulator:
    '''Efficient coverage accumulator using numpy arrays'''
    def __init__(self, length: int, stranded: bool = False):
        self.stranded = stranded
        if stranded:
            self.forward = np.zeros(length, dtype=np.int32)
            self.reverse = np.zeros(length, dtype=np.int32)
        else:
            self.combined = np.zeros(length, dtype=np.int32)
    
    def add_coverage(self, position: int, count: int, is_reverse: bool = False):
        '''Add coverage at a specific position'''
        if self.stranded:
            if is_reverse:
                self.reverse[position] += count
            else:
                self.forward[position] += count
        else:
            self.combined[position] += count

    def get_coverage_chunks(self, chunk_size: int = None):
        '''Generator yielding non-zero coverage in chunks'''
        if chunk_size is None:
            chunk_size = write_buffer  # Use global write_buffer value
        '''Generator yielding non-zero coverage in chunks'''
        if self.stranded:
            for strand, array in [('forward', self.forward), ('reverse', self.reverse)]:
                nonzero = np.nonzero(array)[0]
                for chunk_start in range(0, len(nonzero), chunk_size):
                    chunk_indices = nonzero[chunk_start:chunk_start + chunk_size]
                    yield strand, chunk_indices, array[chunk_indices]
        else:
            nonzero = np.nonzero(self.combined)[0]
            for chunk_start in range(0, len(nonzero), chunk_size):
                chunk_indices = nonzero[chunk_start:chunk_start + chunk_size]
                yield 'combined', chunk_indices, self.combined[chunk_indices]

def process_reads_chunk(chunk_reads, length_offsets: Dict[int, int], accumulator: CoverageAccumulator):
    '''Process a chunk of reads efficiently'''
    # Pre-calculate valid read lengths
    valid_lengths = set(length_offsets.keys())
    
    for read in chunk_reads:
        read_length = read.qlen
        
        # Skip invalid reads
        if read_length < 25 or read_length not in valid_lengths:
            continue
            
        # Get positions and offset
        positions = np.array(read.positions)
        offset = length_offsets[read_length]
        
        # Calculate read count
        read_count = int(read.qname.split("_x")[1]) if "_x" in read.qname else 1
        
        # Calculate A-site
        if not read.is_reverse:
            Asite = positions[offset]
        else:
            Asite = positions[-1 - offset]
            
        accumulator.add_coverage(Asite, read_count, read.is_reverse)

def write_coverage_chunk(chrom: str, indices: np.ndarray, values: np.ndarray, outfile: str):
    '''Write a chunk of coverage data efficiently'''
    with open(outfile, 'a') as f:
        for pos, count in zip(indices, values):
            f.write(f"{chrom}\t{pos}\t{pos+1}\t{count}\n")

def process_bam_file(filepath: str, length_offsets: Dict[int, int], output_prefix: str, stranded: bool = False):
    '''Process BAM file with optimized data structures and I/O'''
    alignments = pysam.AlignmentFile(filepath, "rb")
    
    # Setup output files
    if stranded:
        outfiles = {
            'forward': f"{output_prefix}.forward.bedgraph",
            'reverse': f"{output_prefix}.reverse.bedgraph"
        }
    else:
        outfiles = {
            'combined': f"{output_prefix}.bedgraph"
        }
    
    # Initialize empty files
    for fname in outfiles.values():
        open(fname, 'w').close()
    
    # Process each chromosome
    for chrom in alignments.references:
        try:
            # Get chromosome length from header
            chrom_length = alignments.get_reference_length(chrom)
            logger.debug(f"Processing {chrom} (length: {chrom_length})")
            
            # Initialize coverage accumulator
            accumulator = CoverageAccumulator(chrom_length, stranded)
            
            # Process reads in chunks
            reads_iterator = alignments.fetch(chrom)
            while True:
                chunk = list(itertools.islice(reads_iterator, chunk_size))
                if not chunk:
                    break
                process_reads_chunk(chunk, length_offsets, accumulator)
            
            # Write accumulated coverage
            for strand, indices, values in accumulator.get_coverage_chunks():
                write_coverage_chunk(chrom, indices, values, outfiles[strand])
            
            # Clear memory
            del accumulator
            
        except ValueError as e:
            logger.warning(f"Error processing {chrom}: {e}")
            continue

def main():
    parser = argparse.ArgumentParser(description="Process BAM file using length-specific offsets")
    parser.add_argument("-b", "--bam", 
                      required=True,
                      help="Path to the BAM file (with index)")
    parser.add_argument("-o", "--offsets",
                      required=True,
                      help="Tab-delimited file with read length specific offsets")
    parser.add_argument("-p", "--prefix",
                      required=True,
                      type=Path,
                      help="Output file prefix")
    parser.add_argument("--stranded",
                      action="store_true",
                      help="Generate strand-specific output files")
    parser.add_argument("--chunk-size",
                      type=int,
                      default=DEFAULT_CHUNK_SIZE,
                      help="Number of reads to process at once")
    parser.add_argument("--write-buffer",
                      type=int,
                      default=DEFAULT_WRITE_BUFFER,
                      help="Number of positions to buffer before writing")
    parser.add_argument("--debug",
                      action="store_true",
                      help="Enable debug logging")
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)
        
    # Update global processing parameters if provided
    global chunk_size, write_buffer
    chunk_size = args.chunk_size
    write_buffer = args.write_buffer
    
    # Create output directory if needed
    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Read offsets
        logger.info("Reading offsets file")
        length_offsets = read_offset_file(args.offsets)
        logger.info(f"Loaded offsets for {len(length_offsets)} read lengths")
        
        # Process BAM and write output(s)
        logger.info(f"Processing BAM file {'with' if args.stranded else 'without'} strand separation")
        process_bam_file(args.bam, length_offsets, args.prefix, args.stranded)
        
        logger.info("Processing completed successfully")
        
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        raise

if __name__ == "__main__":
    main()