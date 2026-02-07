# Run Statistics

## Overview

The `RUN_STATISTICS` process generates comprehensive gene annotation statistics for Ensembl species homepages, computing counts and metrics for genes, transcripts, proteins, and other genomic features.

## Process Details

- **Label**: `fetch_file`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Publish Directory**: `${params.outdir}/${meta.gca}`
- **Max Forks**: 20 (limits parallel execution)

## Inputs

| Name | Type | Description |
|------|------|-------------|
| meta | val | Metadata map containing genome assembly and database connection details |

### Required Metadata Fields

- `gca`: Genome assembly accession
- `dbname`: Database name
- `production_name`: Species production name

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| statistics_output | tuple | Metadata and generated SQL files from `core_statistics/` directory |
| versions_file | path | versions.yml file tracking Perl version |

## Parameters

### Required
- `params.outdir`: Output directory for results
- `params.enscode`: Path to Ensembl code repository
- `params.host`: Database host
- `params.port`: Database port

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process:
1. Executes `generate_species_homepage_stats.pl` Perl script
2. Connects to the specified Ensembl core database
3. Computes comprehensive annotation statistics including:
   - Gene counts (by biotype)
   - Transcript counts
   - Protein-coding gene metrics
   - Alternative splicing statistics
   - Non-coding RNA counts
   - Other genomic feature counts
4. Generates SQL files for inserting statistics into the database
5. Outputs results to `core_statistics/` subdirectory
6. Captures Perl version information

## Dependencies

- Perl 5
- `generate_species_homepage_stats.pl` script (from `ensembl-genes/src/perl/ensembl/genes/`)
- Ensembl Perl API
- Database read access

## Generated Statistics

The SQL files typically include statistics for:

### Gene Metrics
- Total gene count
- Protein-coding genes
- Pseudogenes
- Non-coding RNA genes (by type: lncRNA, miRNA, snRNA, etc.)

### Transcript Metrics
- Total transcript count
- Protein-coding transcripts
- Alternative splicing frequency

### Protein Metrics
- Protein-coding sequences
- Average protein length
- Protein domains

### Other Features
- Alternative sequences
- Gene predictions
- Imported annotations

## Notes

- Results are published to a genome-specific subdirectory
- Maximum of 20 concurrent processes to prevent database overload
- The process includes a configurable sleep delay after completion to handle file system latency
- Generated SQL files are designed for the Ensembl core database schema
- Statistics are used to populate species homepage displays in Ensembl browser
- Can be executed against production or staging databases