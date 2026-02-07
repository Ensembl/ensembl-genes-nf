# CLEANING Module

## Overview

The `CLEANING` module is a utility process that removes temporary files and directories created during pipeline execution. It's typically used at the end of workflows to clean up intermediate data and free disk space.

**Module Location**: `pipelines/statistics/modules/cleaning.nf`

## Functionality

This module performs a simple but important function:

- **Removes output directories** for specific taxon/run combinations
- **Frees disk space** after analysis completion
- **Cleans temporary files** that are no longer needed

The module is optional and only used when disk space management is critical.

## Inputs

### Channel Input

```groovy
tuple val(taxon_id), val(run_accession)
```

**Input Parameters**:
- `taxon_id`: NCBI taxonomy identifier (e.g., "9606")
- `run_accession`: Run identifier or accession (e.g., "SRR123456")

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `params.outDir` | String | `./results` | Base output directory to clean |

## Outputs

### Channel Outputs

**None** - This module does not emit any channels.

### Side Effects

**Deleted Directory**: `${params.outDir}/${taxon_id}/${run_accession}/`

All contents of this directory are permanently removed.

## Process Configuration

### Directives

```groovy
scratch false           // Don't use scratch directory
label 'default'         // Use default resource allocation
tag "cleaning"          // Tag for identification
```

### Resource Allocation

From `nextflow.config` (`default` label):

- **CPUs**: 1
- **Memory**: 2 GB
- **Time**: 1 hour
- **Queue**: Standard

### Container

```
ensemblorg/ensembl-genes-metadata:latest
```

**Required Tools**: Standard Unix utilities (`rm`, `bash`)

## Implementation Details

### Cleanup Logic

The module executes a simple `rm -rf` command:

```bash
rm -rf ${params.outDir}/${taxon_id}/${run_accession}
```

**WARNING**: This command permanently deletes all files and subdirectories without confirmation.

### Safety Considerations

The module uses `joinPath()` function in the original implementation, which may cause issues. The actual execution uses string interpolation:

```bash
rm -rf ${params.outDir}/${taxon_id}/${run_accession}
```

**Current Issue**: The original code has a bug:
```groovy
// INCORRECT - joinPath not in script block context
script:
"""
rm -rf joinPath(params.outDir, "${taxon_id}", "${run_accession}")
"""
```

**Corrected version** should be:
```groovy
script:
def targetDir = file(params.outDir).resolve(taxon_id).resolve(run_accession)
"""
rm -rf ${targetDir}
"""
```

## Usage Example

### In a Workflow

```groovy
include { CLEANING } from '../modules/cleaning.nf'

workflow {
    // Create cleanup channel
    def cleanup_ch = channel.of(
        ['9606', 'human_analysis_20240101'],
        ['10090', 'mouse_analysis_20240101']
    )
    
    // Run cleanup
    CLEANING(cleanup_ch)
}
```

### After Analysis Completion

```groovy
workflow {
    // Run analysis
    def results = SOME_ANALYSIS(input_data)
    
    // Extract taxon and run info for cleanup
    def cleanup_ch = results.map { meta, files ->
        tuple(meta.taxon_id, meta.run_accession)
    }
    
    // Clean up when done
    CLEANING(cleanup_ch)
}
```

## Error Handling

### Common Errors

#### 1. Directory Not Found

**Error Message**:
```
rm: cannot remove 'results/9606/run123': No such file or directory
```

**Behavior**: Process completes with exit code 1 (failure)

**Solution**: 
- Check if directory exists before cleanup
- Use conditional cleanup logic

#### 2. Permission Denied

**Error Message**:
```
rm: cannot remove 'results/9606/run123': Permission denied
```

**Solution**:
- Verify file ownership
- Check directory permissions
- Run pipeline with appropriate user privileges

#### 3. Directory in Use

**Error Message**:
```
rm: cannot remove 'results/9606/run123': Device or resource busy
```

**Solution**:
- Wait for other processes to complete
- Check for open file handles
- Unmount any mounted filesystems

## Best Practices

### 1. Conditional Cleanup

Only clean directories that exist:

```groovy
process CLEANING {
    input:
    val taxon_id
    val run_accession

    script:
    def targetDir = "${params.outDir}/${taxon_id}/${run_accession}"
    """
    if [ -d "${targetDir}" ]; then
        echo "Cleaning ${targetDir}..."
        rm -rf "${targetDir}"
    else
        echo "Directory ${targetDir} not found, skipping cleanup"
    fi
    """
}
```

### 2. Safety Checks

Add safety checks before deletion:

```groovy
script:
"""
# Verify path is under outDir (prevent accidental deletion of wrong directories)
TARGET="${params.outDir}/${taxon_id}/${run_accession}"

# Safety check: ensure we're in the output directory
if [[ "\${TARGET}" == ${params.outDir}/* ]]; then
    rm -rf "\${TARGET}"
else
    echo "ERROR: Target path outside output directory"
    exit 1
fi
"""
```

### 3. Dry Run Mode

Add dry-run capability for testing:

```groovy
script:
def dryRun = params.dryRun ? "echo '[DRY RUN]'" : ""
"""
${dryRun} rm -rf ${params.outDir}/${taxon_id}/${run_accession}
"""
```

Usage:
```bash
nextflow run main.nf --dryRun true
```

### 4. Cleanup Logging

Log what's being deleted:

```groovy
script:
"""
TARGET="${params.outDir}/${taxon_id}/${run_accession}"
if [ -d "\${TARGET}" ]; then
    SIZE=\$(du -sh "\${TARGET}" | cut -f1)
    echo "Cleaning \${TARGET} (size: \${SIZE})"
    rm -rf "\${TARGET}"
    echo "Cleanup complete"
else
    echo "Directory \${TARGET} not found"
fi
"""
```

## Performance Considerations

### Execution Time

Deletion time depends on:
- Number of files
- File sizes
- Filesystem type

**Typical deletion times**:
- Small directory (<1 GB, <1000 files): 1-5 seconds
- Medium directory (1-10 GB, 1000-10000 files): 5-30 seconds
- Large directory (>10 GB, >10000 files): 30-300 seconds

### Parallelization

Multiple cleanup operations can run in parallel:

```groovy
// Clean 10 directories simultaneously
channel.of(
    ['9606', 'run1'],
    ['10090', 'run2'],
    // ... 8 more
)
| CLEANING  // Runs 10 instances in parallel
```

**Disk I/O consideration**: Too many parallel deletions can saturate disk I/O.

**Recommended concurrency**: 5-10 parallel deletions

## Advanced Usage

### Selective Cleanup

Clean only specific file types:

```groovy
process CLEANING_SELECTIVE {
    input:
    val taxon_id
    val run_accession
    
    script:
    """
    # Remove only temporary files, keep results
    find ${params.outDir}/${taxon_id}/${run_accession} -name "*.tmp" -delete
    find ${params.outDir}/${taxon_id}/${run_accession} -name "work*" -type d -exec rm -rf {} +
    """
}
```

### Age-Based Cleanup

Remove old directories:

```groovy
process CLEANUP_OLD {
    script:
    """
    # Remove directories older than 30 days
    find ${params.outDir} -maxdepth 2 -type d -mtime +30 -exec rm -rf {} +
    """
}
```

### Archive Before Cleanup

Archive data before deletion:

```groovy
process ARCHIVE_AND_CLEAN {
    input:
    val taxon_id
    val run_accession
    
    output:
    path("*.tar.gz")
    
    script:
    """
    # Archive directory
    tar -czf ${taxon_id}_${run_accession}.tar.gz ${params.outDir}/${taxon_id}/${run_accession}
    
    # Remove original
    rm -rf ${params.outDir}/${taxon_id}/${run_accession}
    """
}
```

## Testing

### Unit Test

Test cleanup for a test directory:

```bash
# Create test directory
mkdir -p results/test_taxon/test_run
echo "test data" > results/test_taxon/test_run/test.txt

# Run cleanup
nextflow run pipelines/statistics/main.nf \
    -entry CLEANING \
    --outDir results \
    -profile docker

# Verify deletion
ls results/test_taxon/test_run 2>/dev/null && echo "FAIL: Directory still exists" || echo "PASS: Directory deleted"
```

## Troubleshooting

### Debug Mode

Check what would be deleted without actually deleting:

```groovy
script:
"""
# List files that would be deleted
echo "Would delete:"
ls -lR ${params.outDir}/${taxon_id}/${run_accession}

# Uncomment to actually delete
# rm -rf ${params.outDir}/${taxon_id}/${run_accession}
"""
```

### Check Disk Space Before/After

```groovy
script:
"""
# Show disk space before
echo "Disk space before cleanup:"
df -h ${params.outDir}

# Remove directory
rm -rf ${params.outDir}/${taxon_id}/${run_accession}

# Show disk space after
echo "Disk space after cleanup:"
df -h ${params.outDir}
"""
```

## Alternatives to This Module

### 1. Nextflow's `cleanup` Parameter

Use Nextflow's built-in cleanup:

```bash
# Clean work directory after successful completion
nextflow run main.nf -with-report -with-dag -cleanup
```

### 2. `publishDir` with `mode: 'move'`

Move files instead of copying, automatically cleaning source:

```groovy
process SOME_PROCESS {
    publishDir "${params.outdir}", mode: 'move'
    // Files are moved, not copied, automatically cleaning work directory
}
```

### 3. Cron Job Cleanup

Schedule periodic cleanup:

```bash
# Cron job: clean old results every week
0 0 * * 0 find /path/to/results -mtime +30 -type d -exec rm -rf {} +
```

## Recommendations

### When to Use This Module

✅ **Use when**:
- You need to clean specific analysis directories
- Running multiple iterations and need to free space between runs
- You want explicit control over what gets deleted
- Part of a production workflow with strict cleanup requirements

❌ **Don't use when**:
- Nextflow's built-in `-cleanup` flag is sufficient
- You want to keep all results for auditing
- Running on temporary filesystems that auto-cleanup
- Unsure about what should be deleted

### Better Alternatives

For most use cases, prefer Nextflow's built-in options:

```bash
# Clean work directory after success
nextflow run main.nf -with-cleanup

# Remove work directory manually after confirmation
nextflow clean -f

# Use tmpfs for work directory (auto-cleanup on reboot)
nextflow run main.nf -work-dir /tmp/nf_work
```

## Related Documentation

- [Nextflow Cleanup Documentation](https://www.nextflow.io/docs/latest/cli.html#clean) - Built-in cleanup options
- [Pipeline Overview](../overview.md) - Pipeline architecture
- [Configuration Guide](../configuration.md) - Output directory configuration

---

**Last Updated**: 2026-02-06  
**Module Version**: 1.0.0  
**Maintained By**: Ensembl Genes Team

**⚠️ WARNING**: This module permanently deletes data. Use with caution and ensure proper backups exist.
