# OMAmer HOG

## Overview

The `OMAMER_HOG` process performs orthology inference using OMAmer (Orthologous MAtrix MER), searching protein sequences against a hierarchical orthologous groups (HOG) database to assign proteins to orthologs.

## Process Details

- **Label**: `omamer`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Store Directory**: `${params.cacheDir}/${meta.gca}/omamer/`
- **Max Forks**: 15 (limits parallel execution)

## Inputs

| Name | Type | Description |
|------|------|-------------|
| meta | val | Metadata map containing genome assembly info |
| translation_file | path | FASTA file containing protein translations to analyze |

### Required Metadata Fields

- `gca`: Genome assembly accession

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| omamer_hog_output | tuple | Metadata and OMAmer output file (`proteins.omamer`) |
| versions_file | path | versions.yml file tracking OMAmer version |

## Parameters

### Required
- `params.cacheDir`: Cache directory for storing OMAmer results
- `params.omamer_database`: Path to OMAmer HOG database

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process:
1. Executes `omamer search` against the specified HOG database
2. Analyzes input protein sequences from the translation file
3. Assigns proteins to hierarchical orthologous groups
4. Outputs results to `proteins.omamer` file
5. Captures OMAmer version from pip package information

## Dependencies

- OMAmer (Python package)
- OMAmer HOG database
- Python 3 with pip

## Storage Strategy

Uses `storeDir` directive to cache results permanently by genome assembly accession, avoiding redundant orthology searches for the same proteins.

## Output Format

The `proteins.omamer` file contains orthology assignments in OMAmer's native format, which can be consumed by downstream tools like OMark for quality assessment.

## Notes

- Results are cached permanently to avoid re-running expensive orthology searches
- Maximum of 15 concurrent processes to manage computational resources
- The process includes a configurable sleep delay after completion to handle file system latency
- OMAmer version is extracted from pip package metadata