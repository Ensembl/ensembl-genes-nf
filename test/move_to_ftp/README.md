# MOVE_TO_FTP Test Pipeline

This test pipeline demonstrates and tests the `MOVE_TO_FTP` module functionality.

## What it does

1. Creates a test file with timestamp and hostname information
2. Transfers the file to an FTP destination using the MOVE_TO_FTP module

## Prerequisites

- Access to the datamover queue (SLURM)
- `become genebuild` command available
- Write access to the FTP destination

## Usage

### Local testing (without SLURM)

```bash
cd test/move_to_ftp
nextflow run main.nf --ftp_destination /path/to/destination
```

### SLURM testing (with datamover queue)

```bash
cd test/move_to_ftp
nextflow run main.nf -profile slurm --ftp_destination /path/to/ftp/destination
```

### Testing with recursive copy

```bash
nextflow run main.nf --ftp_destination /path/to/destination
```

Then edit `nextflow.config` to uncomment `ext.args = '-r'` for recursive copy.

## Parameters

- `--ftp_destination`: Path where files should be copied (required, default: `/tmp/ftp_test_destination`)
- `--outdir`: Output directory for results (default: `results`)

## Expected Output

After successful execution:
- Test file should be present at the FTP destination
- Versions file in the results directory
- Nextflow work directory with execution logs

## Testing the module

Before running on production data:

1. Test locally first to validate the workflow logic:
   ```bash
   nextflow run main.nf
   ```

2. Test with `-stub` to validate without actual execution:
   ```bash
   nextflow run main.nf -stub
   ```

3. Test on SLURM with real datamover queue:
   ```bash
   nextflow run main.nf -profile slurm --ftp_destination /your/ftp/path
   ```
