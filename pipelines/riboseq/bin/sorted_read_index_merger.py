#!/usr/bin/env python3

import gzip
import polars as pl
from pathlib import Path
import time
import bisect
from typing import Dict, List, Tuple, Set, Optional
import json
import tempfile
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

class BatchIndex:
    """Represents a processed batch of samples with their unique reads"""
    
    def __init__(self, batch_id: int, temp_dir: Path):
        self.batch_id = batch_id
        self.temp_dir = temp_dir
        self.unique_reads = []  # Sorted list of sequences
        self.read_to_id = {}    # sequence -> local_id mapping
        self.count_matrix = {}  # local_id -> {sample: count}
        self.samples = []
        self.next_id = 0
        
    def save_to_disk(self):
        """Save batch index to disk and free memory"""
        batch_file = self.temp_dir / f"batch_{self.batch_id}.json"
        
        data = {
            'batch_id': self.batch_id,
            'unique_reads': self.unique_reads,
            'read_to_id': self.read_to_id,
            'count_matrix': self.count_matrix,
            'samples': self.samples,
            'next_id': self.next_id
        }
        
        with open(batch_file, 'w') as f:
            json.dump(data, f)
        
        # Clear memory
        self.unique_reads = []
        self.read_to_id = {}
        self.count_matrix = {}
        
        return batch_file
    
    @classmethod
    def load_from_disk(cls, batch_file: Path):
        """Load batch index from disk"""
        with open(batch_file, 'r') as f:
            data = json.load(f)
        
        batch = cls(data['batch_id'], batch_file.parent)
        batch.unique_reads = data['unique_reads']
        batch.read_to_id = data['read_to_id']
        batch.count_matrix = {int(k): v for k, v in data['count_matrix'].items()}  # JSON converts int keys to str
        batch.samples = data['samples']
        batch.next_id = data['next_id']
        
        return batch

def process_single_file(fasta_gz_path: Path) -> Tuple[str, List[Tuple[str, int]], int]:
    """Extract reads and counts from one file, return sorted"""
    reads_counts = {}
    sample_name = fasta_gz_path.stem.split('_')[0]
    
    try:
        with gzip.open(fasta_gz_path, 'rt') as f:
            while True:
                header_line = f.readline()
                if not header_line:
                    break
                if header_line.startswith('>'):
                    count = int(header_line.strip().split('_x')[1])
                    sequence = f.readline().strip()
                    if sequence:
                        reads_counts[sequence] = count
    except Exception as e:
        print(f"Error processing {fasta_gz_path}: {e}")
        return sample_name, [], 0
    
    # Return sorted list of (sequence, count) tuples
    sorted_reads = sorted(reads_counts.items())
    total_reads = sum(count for _, count in sorted_reads)
    
    return sample_name, sorted_reads, total_reads

def process_batch_worker(args) -> Path:
    """Worker function to process a batch of files"""
    batch_id, file_paths, temp_dir = args
    
    print(f"Processing batch {batch_id} with {len(file_paths)} files...")
    
    batch = BatchIndex(batch_id, temp_dir)
    total_reads_in_batch = 0
    
    for i, file_path in enumerate(file_paths, 1):
        sample_name, sorted_reads, total_reads = process_single_file(file_path)
        
        if not sorted_reads:
            continue
        
        total_reads_in_batch += total_reads
        batch.samples.append(sample_name)
        
        # Sequential merge within batch (same logic as original)
        if not batch.unique_reads:
            # First sample in batch
            for sequence, count in sorted_reads:
                local_id = batch.next_id
                batch.next_id += 1
                
                batch.unique_reads.append(sequence)
                batch.read_to_id[sequence] = local_id
                batch.count_matrix[local_id] = {sample_name: count}
        else:
            # Merge with existing batch sequences
            merge_sorted_reads_into_batch(batch, sorted_reads, sample_name)
    
    compression_ratio = total_reads_in_batch / len(batch.unique_reads) if batch.unique_reads else 0
    print(f"Batch {batch_id} complete: {len(batch.unique_reads):,} unique reads, {compression_ratio:.1f}x compression")
    
    # Save to disk and return path
    return batch.save_to_disk()

def merge_sorted_reads_into_batch(batch: BatchIndex, new_reads: List[Tuple[str, int]], sample_name: str):
    """Merge new sorted reads into existing batch (same logic as original merge)"""
    
    i = 0  # pointer for new_reads
    j = 0  # pointer for existing unique_reads
    new_unique_reads = []
    
    while i < len(new_reads) and j < len(batch.unique_reads):
        new_seq, new_count = new_reads[i]
        existing_seq = batch.unique_reads[j]
        
        if new_seq < existing_seq:
            # New sequence - add to unique set
            local_id = batch.next_id
            batch.next_id += 1
            
            new_unique_reads.append((j, new_seq))  # Insert position, sequence
            batch.read_to_id[new_seq] = local_id
            batch.count_matrix[local_id] = {sample_name: new_count}
            i += 1
            
        elif new_seq > existing_seq:
            # Existing sequence, new sample doesn't have it
            j += 1
            
        else:  # new_seq == existing_seq
            # Sequence exists, add count for this sample
            local_id = batch.read_to_id[existing_seq]
            batch.count_matrix[local_id][sample_name] = new_count
            i += 1
            j += 1
    
    # Add any remaining new sequences
    while i < len(new_reads):
        new_seq, new_count = new_reads[i]
        local_id = batch.next_id
        batch.next_id += 1
        
        new_unique_reads.append((len(batch.unique_reads), new_seq))
        batch.read_to_id[new_seq] = local_id
        batch.count_matrix[local_id] = {sample_name: new_count}
        i += 1
    
    # Insert new unique reads into sorted list (from right to left)
    for insert_pos, sequence in reversed(new_unique_reads):
        batch.unique_reads.insert(insert_pos, sequence)

def merge_two_batches(batch1: BatchIndex, batch2: BatchIndex, output_batch_id: int, temp_dir: Path) -> BatchIndex:
    """Merge two batch indices using two-pointer merge"""
    
    merged = BatchIndex(output_batch_id, temp_dir)
    merged.samples = batch1.samples + batch2.samples
    
    i = j = 0  # Pointers for batch1 and batch2 unique_reads
    
    while i < len(batch1.unique_reads) and j < len(batch2.unique_reads):
        seq1 = batch1.unique_reads[i]
        seq2 = batch2.unique_reads[j]
        
        if seq1 < seq2:
            # Take from batch1
            merged.unique_reads.append(seq1)
            merged.read_to_id[seq1] = merged.next_id
            
            # Copy counts from batch1
            batch1_id = batch1.read_to_id[seq1]
            merged.count_matrix[merged.next_id] = batch1.count_matrix[batch1_id].copy()
            
            merged.next_id += 1
            i += 1
            
        elif seq1 > seq2:
            # Take from batch2
            merged.unique_reads.append(seq2)
            merged.read_to_id[seq2] = merged.next_id
            
            # Copy counts from batch2
            batch2_id = batch2.read_to_id[seq2]
            merged.count_matrix[merged.next_id] = batch2.count_matrix[batch2_id].copy()
            
            merged.next_id += 1
            j += 1
            
        else:  # seq1 == seq2
            # Merge counts from both batches
            merged.unique_reads.append(seq1)
            merged.read_to_id[seq1] = merged.next_id
            
            # Combine counts
            batch1_id = batch1.read_to_id[seq1]
            batch2_id = batch2.read_to_id[seq2]
            
            merged_counts = batch1.count_matrix[batch1_id].copy()
            merged_counts.update(batch2.count_matrix[batch2_id])
            merged.count_matrix[merged.next_id] = merged_counts
            
            merged.next_id += 1
            i += 1
            j += 1
    
    # Add remaining sequences from batch1
    while i < len(batch1.unique_reads):
        seq1 = batch1.unique_reads[i]
        merged.unique_reads.append(seq1)
        merged.read_to_id[seq1] = merged.next_id
        
        batch1_id = batch1.read_to_id[seq1]
        merged.count_matrix[merged.next_id] = batch1.count_matrix[batch1_id].copy()
        
        merged.next_id += 1
        i += 1
    
    # Add remaining sequences from batch2
    while j < len(batch2.unique_reads):
        seq2 = batch2.unique_reads[j]
        merged.unique_reads.append(seq2)
        merged.read_to_id[seq2] = merged.next_id
        
        batch2_id = batch2.read_to_id[seq2]
        merged.count_matrix[merged.next_id] = batch2.count_matrix[batch2_id].copy()
        
        merged.next_id += 1
        j += 1
    
    return merged

class HybridUniqueReadsMerger:
    """Hybrid approach: batch processing + hierarchical merge"""

    def __init__(self, output_dir: Path, batch_size: int = 20, n_workers: int = None,
                 previous_index: Path = None, output_suffix: str = ''):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        self.batch_size = batch_size
        self.n_workers = n_workers or min(mp.cpu_count(), 4)  # Conservative default
        self.previous_index = Path(previous_index) if previous_index else None
        self.output_suffix = output_suffix

        # Create temporary directory for batch files
        self.temp_dir = self.output_dir / "temp_batches"
        self.temp_dir.mkdir(exist_ok=True)

        # Stats
        self.total_reads_processed = 0
        self.total_samples_processed = 0
        
    def chunk_files(self, file_paths: List[Path]) -> List[List[Path]]:
        """Split file paths into batches"""
        return [file_paths[i:i + self.batch_size] 
                for i in range(0, len(file_paths), self.batch_size)]
    
    def process_batches_parallel(self, file_paths: List[Path]) -> List[Path]:
        """Process batches in parallel"""
        
        batches = self.chunk_files(file_paths)
        batch_files = []
        
        print(f"Processing {len(batches)} batches with {self.n_workers} workers...")
        
        # Prepare arguments for workers
        worker_args = [(i, batch, self.temp_dir) for i, batch in enumerate(batches)]
        
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit all batch jobs
            future_to_batch = {executor.submit(process_batch_worker, args): args[0] 
                             for args in worker_args}
            
            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch_id = future_to_batch[future]
                try:
                    batch_file = future.result()
                    batch_files.append(batch_file)
                except Exception as e:
                    print(f"Batch {batch_id} failed: {e}")
        
        batch_files.sort()  # Ensure consistent ordering
        return batch_files
    
    def tournament_merge_batches(self, batch_files: List[Path]) -> BatchIndex:
        """Hierarchically merge batch files using tournament approach"""
        
        if len(batch_files) == 1:
            return BatchIndex.load_from_disk(batch_files[0])
        
        print(f"Tournament merge of {len(batch_files)} batches...")
        
        current_batches = [BatchIndex.load_from_disk(bf) for bf in batch_files]
        merge_round = 0
        
        while len(current_batches) > 1:
            merge_round += 1
            next_batches = []
            
            print(f"  Merge round {merge_round}: {len(current_batches)} -> {(len(current_batches) + 1) // 2}")
            
            # Pair up batches and merge
            for i in range(0, len(current_batches), 2):
                if i + 1 < len(current_batches):
                    # Merge pair
                    merged = merge_two_batches(
                        current_batches[i], 
                        current_batches[i + 1],
                        f"merged_r{merge_round}_p{i//2}",
                        self.temp_dir
                    )
                    next_batches.append(merged)
                else:
                    # Odd one out, carry forward
                    next_batches.append(current_batches[i])
            
            current_batches = next_batches
        
        return current_batches[0]
    
    def save_final_outputs(self, final_batch: BatchIndex):
        """Save the final merged results"""

        # 1. Save unique reads FASTA
        fasta_file = self.output_dir / f"unique_reads{self.output_suffix}.fa"
        print(f"Saving {len(final_batch.unique_reads):,} unique reads to FASTA...")

        with open(fasta_file, 'w') as f:
            for sequence in final_batch.unique_reads:
                unique_id = final_batch.read_to_id[sequence]
                f.write(f">unique_read_{unique_id}\n{sequence}\n")

        # 2. Save count matrix
        matrix_file = self.output_dir / f"count_matrix{self.output_suffix}.parquet"
        print(f"Saving count matrix...")

        matrix_data = []
        for unique_id in range(final_batch.next_id):
            row = {'unique_id': unique_id}

            # Add counts for each sample
            for sample in final_batch.samples:
                count = final_batch.count_matrix[unique_id].get(sample, 0)
                row[sample] = count

            matrix_data.append(row)

        df = pl.DataFrame(matrix_data)
        df.write_parquet(matrix_file)

        # 3. Save read mapping
        mapping_file = self.output_dir / f"read_mapping{self.output_suffix}.json"
        with open(mapping_file, 'w') as f:
            json.dump(final_batch.read_to_id, f)

        print(f"Matrix shape: {len(matrix_data):,} unique reads × {len(final_batch.samples)} samples")

        return fasta_file, matrix_file, mapping_file
    
    def calculate_and_save_stats(self, final_batch: BatchIndex, processing_time: float):
        """Calculate and save processing statistics"""
        
        # Calculate totals
        total_reads = sum(
            sum(counts.values()) 
            for counts in final_batch.count_matrix.values()
        )
        
        final_compression = total_reads / len(final_batch.unique_reads)
        storage_reduction = (1 - len(final_batch.unique_reads) / total_reads) * 100
        
        summary = {
            'processing_summary': {
                'total_files_processed': len(final_batch.samples),
                'total_reads_processed': total_reads,
                'total_unique_reads': len(final_batch.unique_reads),
                'compression_ratio': final_compression,
                'storage_reduction_percent': storage_reduction,
                'processing_time_seconds': processing_time,
                'batch_size': self.batch_size,
                'n_workers': self.n_workers
            },
            'samples_processed': final_batch.samples
        }
        
        summary_file = self.output_dir / f"processing_summary{self.output_suffix}.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"\n=== FINAL SUMMARY ===")
        print(f"Files processed: {len(final_batch.samples)}")
        print(f"Total reads: {total_reads:,}")
        print(f"Unique reads: {len(final_batch.unique_reads):,}")
        print(f"Compression ratio: {final_compression:.1f}x")
        print(f"Storage reduction: {storage_reduction:.1f}%")
        print(f"Processing time: {processing_time:.1f} seconds")
        
        return summary_file
    
    def load_previous_index(self) -> Optional[BatchIndex]:
        """Load previous index if it exists for progressive mode"""
        if not self.previous_index or not self.previous_index.exists():
            return None

        print(f"Loading previous index from {self.previous_index}...")

        # Load previous data
        fasta_file = self.previous_index / "unique_reads.fa"
        matrix_file = self.previous_index / "count_matrix.parquet"
        mapping_file = self.previous_index / "read_mapping.json"

        if not all(f.exists() for f in [fasta_file, matrix_file, mapping_file]):
            print(f"Warning: Previous index incomplete, starting fresh")
            return None

        # Create a batch index from previous data
        prev_batch = BatchIndex(batch_id=-1, temp_dir=self.temp_dir)

        # Load unique reads from FASTA
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                sequence = line.strip()
                if sequence:
                    prev_batch.unique_reads.append(sequence)

        # Load read mapping
        with open(mapping_file, 'r') as f:
            prev_batch.read_to_id = json.load(f)

        # Load count matrix
        df = pl.read_parquet(matrix_file)
        prev_batch.samples = [col for col in df.columns if col != 'unique_id']

        # Reconstruct count matrix
        for row in df.iter_rows(named=True):
            unique_id = row['unique_id']
            prev_batch.count_matrix[unique_id] = {
                sample: row[sample] for sample in prev_batch.samples if row.get(sample, 0) > 0
            }

        prev_batch.next_id = len(prev_batch.unique_reads)

        print(f"Loaded previous index: {len(prev_batch.unique_reads):,} unique reads, {len(prev_batch.samples)} samples")

        return prev_batch

    def cleanup_temp_files(self):
        """Clean up temporary batch files"""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def run_pipeline(self, file_paths: List[Path]) -> Dict:
        """Run the complete hybrid pipeline"""
        start_time = time.time()

        print(f"=== HYBRID UNIQUE READS MERGER ===")
        print(f"Files to process: {len(file_paths)}")
        print(f"Batch size: {self.batch_size}")
        print(f"Workers: {self.n_workers}")
        print(f"Output directory: {self.output_dir}")
        print(f"Progressive mode: {self.previous_index is not None}")

        try:
            # Phase 0: Load previous index if in progressive mode
            prev_batch = self.load_previous_index()

            # Phase 1: Process batches in parallel
            batch_files = self.process_batches_parallel(file_paths)

            if not batch_files:
                raise RuntimeError("No batches processed successfully")

            # Phase 2: Tournament merge of batches
            current_batch = self.tournament_merge_batches(batch_files)

            # Phase 2.5: Merge with previous index if available
            if prev_batch:
                print(f"\nMerging with previous index...")
                final_batch = merge_two_batches(prev_batch, current_batch, "final", self.temp_dir)
                print(f"Final merged index: {len(final_batch.unique_reads):,} unique reads, {len(final_batch.samples)} samples")
            else:
                final_batch = current_batch

            # Phase 3: Save final outputs
            fasta_file, matrix_file, mapping_file = self.save_final_outputs(final_batch)
            
            # Phase 4: Save summary
            processing_time = time.time() - start_time
            summary_file = self.calculate_and_save_stats(final_batch, processing_time)
            
            print(f"\n=== OUTPUTS CREATED ===")
            print(f"1. Unique reads FASTA: {fasta_file}")
            print(f"2. Count matrix: {matrix_file}")
            print(f"3. Read mapping: {mapping_file}")
            print(f"4. Processing summary: {summary_file}")
            
            return {
                'fasta_file': fasta_file,
                'matrix_file': matrix_file,
                'mapping_file': mapping_file,
                'summary_file': summary_file
            }
        
        finally:
            # Always cleanup temp files
            self.cleanup_temp_files()

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Merge collapsed FASTA files to create unique read index and count matrix'
    )
    parser.add_argument(
        '--file-list',
        type=str,
        help='Path to file containing list of FASTA files to process (one per line)'
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        help='Directory containing FASTA files (alternative to --file-list)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=20,
        help='Number of files to process per batch (default: 20)'
    )
    parser.add_argument(
        '--n-workers',
        type=int,
        default=None,
        help='Number of parallel workers (default: min(cpu_count, 4))'
    )
    parser.add_argument(
        '--previous-index',
        type=str,
        help='Path to previous index directory for progressive mode (contains unique_reads.fa, count_matrix.parquet, read_mapping.json)'
    )
    parser.add_argument(
        '--output-suffix',
        type=str,
        default='',
        help='Suffix to append to output filenames (e.g., "_run_1")'
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)

    # Get file paths
    file_paths = []

    if args.file_list:
        # Read from file list
        files_list_path = Path(args.file_list)
        if not files_list_path.exists():
            print(f"File list not found: {files_list_path}")
            return

        with open(files_list_path, 'r') as f:
            file_names = [line.strip() for line in f.readlines() if line.strip()]

        # Files can be absolute paths or relative to current directory
        for name in file_names:
            full_path = Path(name)
            if full_path.exists() and full_path.is_file():
                file_paths.append(full_path)

    elif args.input_dir:
        # Scan directory for .fa.gz files
        input_dir = Path(args.input_dir)
        if not input_dir.exists():
            print(f"Input directory not found: {input_dir}")
            return

        file_paths = list(input_dir.glob("*.fa.gz"))

    else:
        print("Error: Either --file-list or --input-dir must be provided")
        return

    print(f"Found {len(file_paths)} valid files to process")

    if not file_paths:
        print("No valid files found!")
        return

    # Initialize hybrid merger
    merger = HybridUniqueReadsMerger(
        output_dir=output_dir,
        batch_size=args.batch_size,
        n_workers=args.n_workers,
        previous_index=Path(args.previous_index) if args.previous_index else None,
        output_suffix=args.output_suffix
    )

    # Run the pipeline
    results = merger.run_pipeline(file_paths)

    print(f"\nHybrid pipeline complete!")
    print(f"Results in: {output_dir}")

if __name__ == "__main__":
    main()