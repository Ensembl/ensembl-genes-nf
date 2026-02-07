# Fetch Proteins

## Overview

The `FETCH_PROTEINS` process retrieves protein translations from an Ensembl core database, generating FASTA files of protein sequences for downstream analysis.

## Process Details

- **Label**: `fetch_file`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Store Directory**: `${params.cacheDir}/${meta.gca}/translations/`
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
- `params.enscode`: Path to Ensembl code repository
- `params.host`: Database host
- `params.port`: Database port
- `params.production_name`: Default production name (can be overridden by metadata)

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process:
1. Executes `dump_translations.pl` Perl script from the Ensembl genes repository
2. Connects to the specified Ensembl core database
3. Extracts all protein translations
4. Outputs a FASTA file named `${meta.production_name}.fa`
5. Generates a versions file tracking the Perl version used

## Dependencies

- Perl 5
- `dump_translations.pl` script (from `ensembl-genes/src/perl/ensembl/genes/`)
- Ensembl Perl API
- Database read access

## Storage Strategy

Uses `storeDir` directive to cache results permanently by genome assembly accession, avoiding redundant database queries for the same genome.

## Notes

- Results are cached permanently to reduce database load
- The process includes a configurable sleep delay after completion to handle file system latency
- Maximum of 20 concurrent processes to prevent database overload
- Output filename is based on the production name for consistency