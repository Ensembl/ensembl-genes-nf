# Statistics Pipeline Modules Documentation

This directory contains comprehensive documentation for all modules in the Ensembl genes statistics pipeline.

## Module Overview

The statistics pipeline consists of 13 modules organized into functional categories:

### Data Retrieval Modules
1. **[fetch-genome](modules/fetch-genome.md)** - Retrieves genome sequences from Ensembl core databases
2. **[fetch-proteins](modules/fetch-proteins.md)** - Extracts protein translations from Ensembl databases

### BUSCO Quality Assessment Modules
3. **[busco-dataset](modules/busco-dataset.md)** - Downloads appropriate BUSCO lineage datasets
4. **[busco-lineage](modules/busco-lineage.md)** - Runs BUSCO assessment on genome sequences or on protein translations
5. **[busco-core-metakeys](modules/busco-core-metakeys.md)** - Patches BUSCO metadata into core databases

### Orthology Analysis Modules
6. **[omamer-hog](modules/omamer-hog.md)** - Performs orthology inference using OMAmer
7. **[omark](modules/omark.md)** - Quality assessment of protein annotations using OMark

### Statistics Generation Modules
8. **[run-statistics](modules/run-statistics.md)** - Generates comprehensive annotation statistics
9. **[run-ensembl-meta](modules/run-ensembl-meta.md)** - Generates core database metadata SQL files


### Database Operations Modules
10. **[populate-db](modules/populate-db.md)** - Executes SQL files to populate databases
11. **[db-metadata](modules/db-metadata.md)** - Manages database metadata and versioning


## Pipeline Flow

The typical execution flow of the statistics pipeline:

```
1. Data Retrieval
   └─> FETCH_GENOME
   └─> FETCH_PROTEINS

2. Quality Assessment (Parallel)
   ├─> BUSCO_DATASET
   │   ├─> BUSCO_LINEAGE
            └─> BUSCO_CORE_METAKEYS
   │   └─> BUSCO_LINEAGE
            └─> BUSCO_CORE_METAKEYS
   │
   └─> OMAMER_HOG
       └─> OMARK

3. Statistics Generation
   ├─> RUN_STATISTICS
   └─> RUN_ENSEMBL_META

4. Database Population
   └─> POPULATE_DB

5. Metadata Management
   └─> DB_METADATA

```

## Module Categories by Function

### Quality Control
- **BUSCO Genome Lineage**: Assesses genome completeness
- **BUSCO Protein Lineage**: Assesses proteome completeness
- **OMark**: Validates annotation consistency using orthology

### Statistics & Metrics
- **Run Statistics**: Computes gene/transcript/protein counts by biotype
- **Run Ensembl Meta**: Generates schema and species metadata

### Database Management
- **BUSCO Core Metakeys**: Inserts BUSCO metrics into database
- **Populate DB**: Executes SQL files for statistics insertion
- **DB Metadata**: Manages database versioning and tracking

### Resource Management
- **Cleaning**: Removes cached intermediate files

## Key Dependencies

### External Tools
- **BUSCO** (v5+): Genome/proteome completeness assessment
- **OMAmer**: Orthology inference
- **OMark**: Annotation quality assessment

### Ensembl Dependencies
- **Ensembl Perl API**: Database access and manipulation
- **Ensembl Python libraries**: Metadata generation
- **Ensembl scripts**: Statistics computation and data extraction

### Databases
- **Ensembl Core Database**: Primary target for statistics
- **BUSCO Lineage Datasets**: Reference for completeness assessment
- **OMAmer HOG Database**: Reference for orthology inference

## Common Parameters

Most modules use these common parameters:

### Database Connection
- `params.host`: Database host
- `params.port`: Database port
- `params.user`: Database username (when required)
- `params.password`: Database password (when required)

### Paths
- `params.outdir`: Output directory for published results
- `params.cacheDir`: Cache directory for intermediate files
- `params.enscode`: Path to Ensembl code repository

### Execution Control
- `params.files_latency`: Delay after file operations (file system sync)
- `maxForks`: Limit parallel process execution

## Caching Strategy

Several modules use `storeDir` for persistent caching:

- **fetch-genome**: Caches genome sequences by GCA
- **fetch-proteins**: Caches protein translations by GCA
- **busco-dataset**: Caches BUSCO datasets by lineage
- **omamer-hog**: Caches orthology assignments by GCA

This strategy reduces redundant computations and database queries when processing the same genomes multiple times.

## Conditional Execution

Some modules execute conditionally based on parameters:

- **BUSCO_CORE_METAKEYS**: `params.apply_busco_metakeys`
- **POPULATE_DB**: `params.apply_ensembl_stats` OR `params.apply_ensembl_beta_metakeys`
- **CLEANING**: `params.clean` AND `params.clean_work_dir`

## Output Structure

Results are typically organized by genome assembly accession:

```
${params.outdir}/
└── ${meta.gca}/
    ├── busco_genome_lineage/
    ├── busco_protein_lineage/
    ├── omark_output/
    ├── core_statistics/
    │   └── *.sql
    └── versions.yml
```

## Metadata Requirements

All modules expect metadata maps with these common fields:
- `gca`: Genome assembly accession (GCA_XXXXXXXXX.X)
- `dbname`: Ensembl core database name
- `production_name`: Species production name
- `species_id`: Species identifier (for database operations)

## Documentation Format

Each module documentation includes:
- **Overview**: Purpose and functionality
- **Process Details**: Labels, tags, directives
- **Inputs**: Expected input channels and metadata
- **Outputs**: Generated output channels
- **Parameters**: Required and optional configuration
- **Script Details**: What the process does
- **Dependencies**: External tools and libraries
- **Notes**: Important considerations and best practices

## Version Tracking

All modules generate `versions.yml` files that track:
- Tool versions (BUSCO, OMAmer, OMark, etc.)
- Language versions (Python, Perl)
- Database client versions (MySQL)

These version files support reproducibility and troubleshooting.

## For More Information

- See individual module documentation files for detailed information
- Refer to the main pipeline documentation for workflow orchestration
- Check Ensembl documentation for database schema details
- Consult tool-specific documentation (BUSCO, OMAmer, OMark) for advanced usage