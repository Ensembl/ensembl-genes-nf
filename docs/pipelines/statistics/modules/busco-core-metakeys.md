# BUSCO Core Metakeys

## Overview

The `BUSCO_CORE_METAKEYS` process patches BUSCO metadata into an Ensembl core database by parsing BUSCO summary files and inserting metakeys directly into the database.

## Process Details

- **Label**: `python`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Cache**: Disabled because the process writes to the database
- **Publish Directory**: `${params.outdir}/${meta.gca}`
- **Conditional Execution**: Only runs when `params.apply_busco_metakeys` is true

## Inputs

| Name | Type | Description |
|------|------|-------------|
| meta | val | Metadata map containing genome assembly info and database details |
| summary_file | path | BUSCO summary file to parse for metakeys |

### Required Metadata Fields

- `gca`: Genome assembly accession
- `dbname`: Database name
- `species_id`: Species identifier in the database

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| metakey_json | tuple val(meta), path | BUSCO metakey JSON file published to `${params.outdir}/${meta.gca}` |
| versions_file | path | Versions YAML file tracking Python version |

## Parameters

### Required
- `params.apply_busco_metakeys`: Boolean flag to enable/disable process execution
- `params.outdir`: Output directory for results
- `params.host`: Database host
- `params.port`: Database port
- `params.user`: Database username
- `params.password`: Database password

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process:
1. Executes `busco_metakeys_patch.py` with database connection parameters
2. Parses the BUSCO summary file
3. Creates a BUSCO metakey JSON file in the task working directory
4. Inserts metakeys into the specified Ensembl core database
5. Runs the query directly against the database (`-run_query true`)
6. Publishes the JSON file to `${params.outdir}/${meta.gca}`
7. Generates a versions file tracking the Python version used

## Dependencies

- Python 3
- `busco_metakeys_patch.py` script (from Ensembl genes repository)
- Database access credentials

## Notes

- The process includes a configurable sleep delay after completion to handle file system latency
- JSON results are published to a genome-specific subdirectory
- Direct database modification requires appropriate write permissions
- Nextflow task caching is disabled for this process because it writes BUSCO metakeys to the database
