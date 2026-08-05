# Fetch Proteins

## Overview

The `FETCH_PROTEINS` process retrieves protein translations from an Ensembl core database, generating FASTA files of protein sequences for downstream analysis.

## Process Details

- **Label**: `fetch_file`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Store Directory**: `${params.cacheDir}/${meta.gca}/fasta/`
- **Max Forks**: 20 (limits parallel execution)

## Inputs

| Name | Type | Description |
|------|------|-------------|
| meta | val | Metadata map containing genome assembly and database connection details |

### Required Metadata Fields

- `gca`: Genome assembly accession
- `dbname`: Database name containing protein translations
- `production_name`: Species production name

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| proteins_file | path | FASTA file containing protein translations |
| versions_file | path | versions.yml file tracking Perl version |

## Parameters

### Required
- `params.cacheDir`: Cache directory for storing protein files
- `params.legacy_enscode`: Temporary path to the legacy Ensembl code repository
- `params.host`: Database host
- `params.port`: Database port
- `params.production_name`: Default production name (can be overridden by metadata)

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process:
1. Uses `storeDir` to reuse cached `translations.fa` and `versions.yml` when they exist
2. Uses `meta.protein_file` directly when a protein FASTA is provided
3. Executes `dump_translations.pl` when no stored output is available
4. Outputs a FASTA file named `translations.fa`
5. Generates a versions file tracking the Perl version used

## Dependencies

- Perl 5
- `dump_translations.pl` script (from `ensembl-genes/src/perl/ensembl/genes/`)
- Ensembl Perl API
- Database read access

## Storage Strategy

Uses the Nextflow `storeDir` directive to cache results permanently by genome assembly accession, avoiding redundant database queries for the same genome.

## Notes

- Results are cached permanently to reduce database load
- The process includes a configurable sleep delay after completion to handle file system latency
- Maximum of 20 concurrent processes to prevent database overload
- Output filename is `translations.fa` for consistency with downstream BUSCO inputs
