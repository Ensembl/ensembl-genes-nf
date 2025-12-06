# Unique Reads Tracking Integration

This document summarizes the integration of unique reads tracking into the riboseq workflow.

## Overview

The riboseq pipeline now automatically tracks all unique read sequences and their occurrence counts across all samples in a run. This creates a comprehensive index that can be used for cross-sample analysis, compression, and progressive database building.

## What Was Added

### 1. New Module: `MERGE_UNIQUE_READS`

**Location:** [pipelines/riboseq/modules/merge_unique_reads.nf](modules/merge_unique_reads.nf)

This Nextflow process:
- Collects all collapsed FASTA files from the run
- Uses the `sorted_read_index_merger.py` script to create unique read index
- Supports two modes: `per_run` and `progressive`
- Outputs: unique_reads.fa, count_matrix.parquet, read_mapping.json, processing_summary.json

### 2. Enhanced Python Script: `sorted_read_index_merger.py`

**Location:** [pipelines/riboseq/bin/sorted_read_index_merger.py](bin/sorted_read_index_merger.py)

Enhancements made:
- **Command-line interface**: Accepts file lists, batch size, worker count
- **Progressive mode**: Can merge with previous index to grow database incrementally
- **Output suffixing**: Supports naming output files with run identifiers
- **Memory efficient**: Processes files in batches with parallel workers
- **Tournament merge**: Hierarchically merges batches for scalability

### 3. Updated POST_PROCESSING Subworkflow

**Location:** [pipelines/riboseq/subworkflows/post_processing.nf](subworkflows/post_processing.nf)

Changes:
- Added `collapsed_fastas` input parameter
- Collects all collapsed FASTA files from samples
- Calls `MERGE_UNIQUE_READS` process
- Exposes unique reads outputs in emit block

### 4. Updated Main Workflow

**Location:** [pipelines/riboseq/main.nf](main.nf)

Changes:
- Passes collapsed FASTA files from `DATA_ACQUISITION` to `POST_PROCESSING`
- Updated comment to reflect new functionality

### 5. Container Support

**Location:** [pipelines/riboseq/containers/unique-reads-merger.Dockerfile](containers/unique-reads-merger.Dockerfile)

Provides:
- Dockerfile for building container with Python 3.10 and polars
- Instructions for building and deploying container
- Conda specification in module for automatic Wave builds

### 6. Documentation

**Location:** [pipelines/riboseq/docs/unique_reads_tracking.md](docs/unique_reads_tracking.md)

Comprehensive documentation covering:
- Overview of unique reads tracking
- Output file descriptions
- Operating modes (per-run vs progressive)
- Container setup options
- Usage examples
- Performance tuning
- Troubleshooting

### 7. Test Script

**Location:** [pipelines/riboseq/test_unique_reads.sh](test_unique_reads.sh)

Provides:
- Automated testing of both per-run and progressive modes
- Verification of output files
- Integration testing

## Output Files

After running the pipeline, you'll find in `results/unique_reads/`:

1. **`unique_reads.fa`** (or `unique_reads_run_<name>.fa`)
   - FASTA file with all unique sequences
   - IDs: `>unique_read_0`, `>unique_read_1`, etc.

2. **`count_matrix.parquet`**
   - Parquet table: `unique_id | sample1 | sample2 | ...`
   - Sparse matrix with counts per sample

3. **`read_mapping.json`**
   - JSON: `{"ATCG...": 0, "GCTA...": 1, ...}`
   - Sequence → unique_id mapping

4. **`processing_summary.json`**
   - Statistics: files processed, compression ratio, timing

## Usage

### Basic Usage (Per-Run Mode)

```bash
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet samples.csv \
  --star_index /path/to/star \
  --gtf /path/to/genes.gtf \
  --fasta /path/to/genome.fa \
  --chrom_sizes_file /path/to/chrom.sizes \
  --unique_reads_mode per_run
```

Each run creates an independent index with suffix `_run_<runname>`.

### Progressive Mode (Growing Database)

```bash
# First batch
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet batch1.csv \
  --unique_reads_mode progressive \
  --outdir ./results_batch1 \
  ...

# Second batch (merges with previous)
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet batch2.csv \
  --unique_reads_mode progressive \
  --unique_reads_previous_index ./results_batch1/unique_reads \
  --outdir ./results_batch2 \
  ...
```

The index grows with each run, maintaining all previous samples.

## Parameters

### `--unique_reads_mode`
- **Default:** `per_run`
- **Options:** `per_run`, `progressive`
- Operating mode for unique reads tracking

### `--unique_reads_previous_index`
- **Required for:** Progressive mode
- Path to directory with previous index files
- Example: `./results/unique_reads`

## Performance

The merger uses:
- **Batch processing**: Processes files in configurable batches
- **Parallel workers**: Multi-core processing within batches
- **Tournament merge**: Hierarchical merging for scalability
- **Memory efficiency**: Streaming and batching to avoid OOM

Default settings (configurable):
- Batch size: 20 files
- Workers: min(CPU cores, 4)

For large datasets, tune in config:

```groovy
process {
  withName: 'MERGE_UNIQUE_READS' {
    cpus = 16
    memory = '64.GB'
    ext {
      batch_size = 30
      n_workers = 16
    }
  }
}
```

## Algorithm

1. **Batch Processing**: Files divided into batches, processed in parallel
2. **Within-batch merge**: Reads sorted and merged sequentially
3. **Tournament merge**: Batches hierarchically merged in pairs
4. **Progressive merge** (optional): Result merged with previous index
5. **Output**: Unified index with all unique reads + count matrix

Complexity: O(N log N) where N = total reads across all samples

## Integration Points

The unique reads tracking integrates at the POST_PROCESSING stage, after:
- ✅ Data acquisition
- ✅ Quality control
- ✅ Alignment
- ✅ RiboMetric/RiboWaltz analysis

This ensures all samples are processed consistently before indexing.

## Compatibility

- **Nextflow:** DSL2
- **Python:** 3.10+
- **Dependencies:** polars 0.20.0
- **Containers:** Docker, Singularity
- **Platforms:** Linux, macOS (with adjustments)

## Future Enhancements

Potential improvements:
- [ ] Support for other input formats (BAM, BED)
- [ ] Built-in visualization of count matrix
- [ ] Integration with downstream analysis tools
- [ ] Compression of count matrix for very large datasets
- [ ] Cloud storage integration (S3, GCS)

## Getting Help

- See [docs/unique_reads_tracking.md](docs/unique_reads_tracking.md) for detailed documentation
- Run `./test_unique_reads.sh` to verify installation
- Check GitHub issues for known problems

## Credits

- **Algorithm**: Hybrid batch processing + hierarchical merge
- **Implementation**: sorted_read_index_merger.py
- **Integration**: EnsEMBL Genebuild team
