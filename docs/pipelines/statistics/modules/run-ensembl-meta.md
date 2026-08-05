# Run Ensembl Meta

## Overview

The `RUN_ENSEMBL_META` process generates SQL files containing Ensembl core database metadata, including schema information, species details, and production metadata required for Ensembl release databases.

## Process Details

- **Label**: `python`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Store Directory**: `${params.cacheDir}/${meta.gca}/core_statistics/metakeys`
- **Publish Directory**: `${params.outdir}/${meta.gca}`

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
| ensembl_meta_output | tuple | Metadata and generated SQL files (`*.sql`) |
| versions_file | path | versions.yml file tracking Python version |

## Parameters

### Required
- `params.cacheDir`: Cache directory for stored process outputs
- `params.outdir`: Output directory for results
- `params.legacy_enscode`: Temporary path to the legacy Ensembl code repository
- `params.host`: Database host
- `params.port`: Database port
- `params.team`: Team identifier for metadata attribution

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process:
1. Executes `core_meta_data.py` Python script from the Ensembl genes repository
2. Connects to the specified core database
3. Generates SQL files with metadata insertions/updates for:
   - Schema version information
   - Species metadata (taxonomy, assembly, etc.)
   - Production database references
   - Team attribution
4. Outputs SQL files to `core_statistics/` subdirectory
5. Copies SQL files into the task working directory for Nextflow output tracking
6. Captures Python version information

## Dependencies

- Python 3
- `core_meta_data.py` script (from `ensembl-genes/src/python/ensembl/genes/metadata/`)
- Ensembl Python libraries
- Database read access

## Generated SQL Content

The SQL files typically include INSERT/UPDATE statements for the `meta` table with keys such as:
- `schema_version`
- `schema_type`
- `assembly.default`
- `species.production_name`
- `species.taxonomy_id`
- `species.scientific_name`

## Notes

- Results are published to a genome-specific subdirectory
- SQL generation is cached with `storeDir`; clear this cache if the source core database or team metadata changes
- The process includes a configurable sleep delay after completion to handle file system latency
- Generated SQL files can be executed using the `POPULATE_DB` process
- The `core_statistics/` directory is created as an output subdirectory
- SQL files are copied into the task directory so they are accessible in the publish directory
