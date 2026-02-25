# Unique Reads Tracking

The riboseq pipeline includes a feature to track unique reads and their occurrences across all samples. This creates a comprehensive index of all unique read sequences and a count matrix showing how many times each unique read appears in each sample.

## Overview

The unique reads tracking system:
- Processes all collapsed FASTA files from your samples
- Identifies every unique read sequence across all samples
- Generates a count matrix showing occurrences of each unique read per sample
- Supports two modes: **per-run** and **progressive**

## Output Files

The process generates four output files in the `unique_reads/` directory:

1. **`unique_reads.fa`** (or `unique_reads_run_<name>.fa` in per-run mode)
   - FASTA file containing all unique read sequences
   - Each sequence has a unique ID: `>unique_read_<id>`

2. **`count_matrix.parquet`**
   - Parquet table with columns: `unique_id`, `<sample1>`, `<sample2>`, ...
   - Each row represents a unique read
   - Columns contain counts for each sample (0 if read not present)

3. **`read_mapping.json`**
   - JSON mapping from read sequence to unique ID
   - Useful for lookups and downstream analysis

4. **`processing_summary.json`**
   - Statistics about the processing:
     - Total files processed
     - Total reads processed
     - Total unique reads found
     - Compression ratio (how many reads collapse to unique sequences)
     - Processing time

## Operating Modes

### Per-Run Mode (Default)

In per-run mode, each pipeline execution creates an independent unique reads index for just the samples in that run.

```bash
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet samples.csv \
  --unique_reads_mode per_run \
  ...
```

Output files will be named with the run identifier:
- `unique_reads_run_<runname>.fa`
- `count_matrix_run_<runname>.parquet`
- etc.

### Progressive Mode

In progressive mode, the index grows across runs by merging new samples with a previous index. This is useful when:
- Processing samples incrementally over time
- Building a comprehensive database across multiple runs
- Maintaining a single unified index

```bash
# First run - creates initial index
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet samples_batch1.csv \
  --unique_reads_mode progressive \
  ...

# Subsequent runs - merges with previous index
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet samples_batch2.csv \
  --unique_reads_mode progressive \
  --unique_reads_previous_index ./results/unique_reads \
  ...
```

The progressive mode will:
1. Load the previous index (unique_reads.fa, count_matrix.parquet, read_mapping.json)
2. Process new samples
3. Merge previous data with new data
4. Save updated index with all samples combined

## Container Setup

The unique reads merger requires Python 3.10 with the `polars` library. You have several options:

### Option 1: Use Conda (Recommended for Development)

The process is configured with conda specifications. If you're using conda/mamba:

```bash
conda create -n unique-reads python=3.10 polars=0.20.0
conda activate unique-reads
```

### Option 2: Build Docker Container

A Dockerfile is provided in `containers/unique-reads-merger.Dockerfile`:

```bash
# Build the container
docker build -t unique-reads-merger:latest \
  -f pipelines/riboseq/containers/unique-reads-merger.Dockerfile \
  pipelines/riboseq/containers/

# Tag for your registry (optional)
docker tag unique-reads-merger:latest your-registry/unique-reads-merger:latest
docker push your-registry/unique-reads-merger:latest
```

Then update the module configuration to use your container:

```groovy
process {
  withName: 'MERGE_UNIQUE_READS' {
    container = 'your-registry/unique-reads-merger:latest'
  }
}
```

### Option 3: Use Wave to Build Container On-Demand

Nextflow can use Wave to build containers automatically from the conda specification:

```bash
nextflow run pipelines/riboseq/main.nf \
  -with-wave \
  -with-conda \
  ...
```

Wave will automatically build a container with the required conda packages.

## Parameters

### `--unique_reads_mode`
- **Type:** String
- **Default:** `per_run`
- **Options:** `per_run`, `progressive`
- **Description:** Operating mode for unique reads tracking

### `--unique_reads_previous_index`
- **Type:** Path
- **Required for:** Progressive mode
- **Description:** Path to directory containing previous index files
- **Example:** `./results/unique_reads`

## Performance Tuning

The unique reads merger uses parallel processing and can be tuned via process configuration:

```groovy
process {
  withName: 'MERGE_UNIQUE_READS' {
    cpus = 16           // More CPUs = more parallel workers
    memory = '64.GB'    // Memory for large datasets
    ext {
      batch_size = 20   // Files per batch (adjust based on memory)
      n_workers = 16    // Parallel workers (usually = cpus)
    }
  }
}
```

### Memory Recommendations

- **Small datasets** (<100 samples, <10M reads each): 16-32 GB
- **Medium datasets** (100-500 samples): 32-64 GB
- **Large datasets** (>500 samples, >50M reads each): 64-128 GB

### Batch Size Guidelines

- Larger batch size = less I/O overhead but more memory per batch
- Smaller batch size = more parallelizable but more disk I/O
- Recommended: 10-30 files per batch

## Example Workflows

### Single Run Analysis

```bash
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet all_samples.csv \
  --star_index /path/to/star/index \
  --gtf /path/to/annotation.gtf \
  --fasta /path/to/genome.fa \
  --chrom_sizes_file /path/to/chrom.sizes \
  --unique_reads_mode per_run
```

### Incremental Database Building

```bash
# Week 1: Process first batch
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet week1_samples.csv \
  --unique_reads_mode progressive \
  --outdir ./results_week1 \
  ...

# Week 2: Add more samples to the index
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet week2_samples.csv \
  --unique_reads_mode progressive \
  --unique_reads_previous_index ./results_week1/unique_reads \
  --outdir ./results_week2 \
  ...

# Week 3: Continue growing the index
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet week3_samples.csv \
  --unique_reads_mode progressive \
  --unique_reads_previous_index ./results_week2/unique_reads \
  --outdir ./results_week3 \
  ...
```

## Analysis Examples

### Loading the Count Matrix

```python
import polars as pl

# Load count matrix
df = pl.read_parquet("results/unique_reads/count_matrix.parquet")

# Show summary
print(f"Unique reads: {len(df)}")
print(f"Samples: {len(df.columns) - 1}")

# Find reads present in all samples
all_samples = [col for col in df.columns if col != 'unique_id']
universal_reads = df.filter(
    pl.all_horizontal([pl.col(s) > 0 for s in all_samples])
)
print(f"Universal reads (in all samples): {len(universal_reads)}")

# Find sample-specific reads
for sample in all_samples:
    specific = df.filter(
        (pl.col(sample) > 0) &
        (pl.sum_horizontal([pl.col(s) for s in all_samples if s != sample]) == 0)
    )
    print(f"{sample} specific reads: {len(specific)}")
```

### Looking Up Read Sequences

```python
import json

# Load read mapping
with open("results/unique_reads/read_mapping.json") as f:
    read_to_id = json.load(f)

# Look up specific sequence
sequence = "ATCGATCGATCG..."
if sequence in read_to_id:
    unique_id = read_to_id[sequence]
    print(f"Sequence maps to unique_read_{unique_id}")

    # Get counts from matrix
    row = df.filter(pl.col("unique_id") == unique_id)
    print(row)
```

## Algorithm Details

The unique reads merger uses a hybrid approach:

1. **Batch Processing**: Files are divided into batches and processed in parallel
2. **Sorted Merge**: Within each batch, reads are kept sorted for efficient merging
3. **Tournament Merge**: Batches are hierarchically merged using a tournament tree
4. **Progressive Integration**: Previous index is merged as a single large batch

This approach provides:
- **Scalability**: Handles hundreds of samples with billions of reads
- **Memory Efficiency**: Processes data in batches, doesn't load everything at once
- **Performance**: Parallel processing of independent batches
- **Correctness**: Guaranteed to find all unique reads across all samples

## Troubleshooting

### Out of Memory Errors

If you encounter OOM errors:
1. Reduce `batch_size` in process configuration
2. Reduce `n_workers` to use fewer parallel processes
3. Increase allocated memory for the process

### Slow Performance

If processing is slow:
1. Increase `n_workers` up to available CPU cores
2. Increase `batch_size` if you have sufficient memory
3. Use local storage instead of network filesystems

### Progressive Mode Not Loading Previous Index

Check that:
1. `--unique_reads_previous_index` points to correct directory
2. Directory contains all three files: `unique_reads.fa`, `count_matrix.parquet`, `read_mapping.json`
3. Files are not corrupted (try loading them manually)

## Technical Notes

- The count matrix uses Parquet format for efficient storage and fast loading
- Unique IDs are assigned sequentially as reads are discovered
- In progressive mode, IDs are preserved from the previous index
- The algorithm guarantees deterministic results (same input = same output)
- Complexity: O(N log N) where N is total number of reads across all samples
