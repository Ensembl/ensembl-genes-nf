# Statistics Pipeline

The **Statistics Pipeline** generates comprehensive quality metrics and statistics for gene annotations and assemblies. This pipeline is essential for validating annotations, assessing completeness, and providing metadata for Ensembl databases.

## Quick Links

- **[Quick Start Guide](quickstart.md)** - Get up and running quickly
- **[Input Specification](input.md)** - Prepare your input data
- **[Output Reference](output.md)** - Understand the results
- **[Troubleshooting](troubleshooting.md)** - Solve common issues

## Workflow Documentation

### Quality Assessment Workflows

| Workflow | Purpose | Key Features |
|----------|---------|--------------|
| **[BUSCO](workflows/busco.md)** | Assess annotation/assembly completeness | Single-copy ortholog presence, protein & genome modes |
| **[OMArk](workflows/omark.md)** | Proteome quality & contamination screening | Consistency checks, lineage validation |
| **[Ensembl Stats](workflows/ensembl-stats.md)** | Generate database statistics | Gene counts, transcript metrics, metakeys |

### Workflow Selection Guide

```
Need to assess...
├─ Assembly quality?
│  └─ Use BUSCO (genome mode)
├─ Annotation completeness?
│  └─ Use BUSCO (protein mode)
├─ Contamination?
│  └─ Use OMArk
├─ Database statistics?
│  └─ Use Ensembl Stats
└─ Complete QC?
   └─ Use all three workflows
```

## Common Use Cases

### 1. Complete Quality Control

Run all workflows for comprehensive assessment:

```bash
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --busco_mode both \
  --run_omark \
  --run_ensembl_stats \
  --host mysql-server.example.com \
  --user_r ensro \
  --enscode /path/to/ENSCODE \
  --outdir qc_results
```

**Provides:**
- ✅ Assembly completeness (BUSCO genome)
- ✅ Annotation completeness (BUSCO protein)
- ✅ Contamination screening (OMArk)
- ✅ Database statistics (Ensembl Stats)

### 2. Pre-Release Validation

Validate before public release:

```bash
nextflow run main.nf \
  --csvFile release_databases.csv \
  --run_busco_core \
  --busco_mode protein \
  --run_omark \
  --run_ensembl_stats \
  --apply_ensembl_stats \
  --host staging-db.example.com \
  --user ensadmin \
  --password ${DB_PASS} \
  --enscode /path/to/ENSCODE \
  --team genebuild \
  --outdir release_validation
```

### 3. NCBI Assembly Assessment

Download and assess NCBI assemblies:

```bash
# Create CSV with NCBI assembly accessions
cat > ncbi_assemblies.csv << EOF
dbname,species_id,taxon_id,assembly_accession,assembly_name,taxon_name
gca_001234567_core,1,9606,GCA_001234567.1,ASM123456v1,homo_sapiens
gca_002345678_core,1,10090,GCA_002345678.1,ASM234567v1,mus_musculus
EOF

nextflow run main.nf \
  --csvFile ncbi_assemblies.csv \
  --run_busco_ncbi \
  --outdir ncbi_assessment
```

### 4. Comparative Analysis

Compare quality across multiple species:

```bash
# Vertebrate comparison
cat > vertebrates.csv << EOF
dbname,species_id,taxon_id
homo_sapiens_core_110_38,1,9606
mus_musculus_core_110_39,1,10090
gallus_gallus_core_110_7,1,9031
danio_rerio_core_110_11,1,7955
EOF

nextflow run main.nf \
  --csvFile vertebrates.csv \
  --run_busco_core \
  --busco_mode protein \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro \
  --outdir vertebrate_comparison
```

## Pipeline Architecture

```
Input CSV
    │
    ├─── BUSCO Analysis
    │    ├─ Protein mode → Annotation completeness
    │    └─ Genome mode → Assembly completeness
    │
    ├─── OMArk Analysis
    │    ├─ Completeness assessment
    │    └─ Contamination detection
    │
    └─── Ensembl Stats
         ├─ Gene/transcript counts
         ├─ Biotype distributions
         └─ Metakey generation
              │
              └─ Apply to database (optional)
```

## Key Features

### Flexible Input Options

- **Core databases**: Connect to existing Ensembl core databases
- **NCBI assemblies**: Automatically download and assess
- **Mixed sources**: Combine different input types

### Multiple Analysis Modes

- **BUSCO**: Protein, genome, or both modes
- **OMArk**: Proteome-based quality with contamination detection
- **Ensembl Stats**: Comprehensive database metrics

### Database Integration

- **Read-only mode**: Generate statistics without modifying databases
- **Apply mode**: Load statistics and metakeys into databases
- **Validation**: Pre-check before applying changes

### Batch Processing

- Process hundreds of genomes in parallel
- Automatic resource management
- Resume capability for interrupted runs

## Output Overview

### Directory Structure

```
results/
├── busco/
│   ├── sample1_busco_short_summary.txt
│   ├── sample1_genome_busco_short_summary.txt
│   └── sample1_busco_full_table.tsv
├── omark/
│   └── sample1_omark_proteins_detailed_summary.txt
└── ensembl_stats/
    └── sample1_statistics.json
```

### Result Interpretation

| Metric | Good Range | Warning Range | Action Needed |
|--------|------------|---------------|---------------|
| **BUSCO Complete** | >95% | 85-95% | <85% |
| **OMArk Consistency** | >98% | 95-98% | <95% |
| **Gene Count** | Expected ±10% | Expected ±20% | Outside ±20% |

## Requirements

### System Requirements

- **Nextflow**: 21.04.0 or higher
- **Java**: 11 or higher
- **Memory**: 32+ GB recommended
- **Storage**: 50+ GB for temporary files

### Software Dependencies

- **BUSCO**: 5.4.0+
- **OMArk**: Latest version
- **Ensembl API**: Release-specific
- **Singularity/Docker**: For containerized workflows

### Database Access

- **MySQL client**: For core database access
- **Read access**: For statistics generation
- **Write access**: For applying metakeys (optional)

## Getting Started

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/
```

### 2. Clone Pipeline

```bash
git clone https://github.com/Ensembl/ensembl-genes.git
cd ensembl-genes/statistics
```

### 3. Prepare Input

Create a CSV file with your targets:

```csv
dbname,species_id,taxon_id
homo_sapiens_core_110_38,1,9606
```

### 4. Run Pipeline

```bash
nextflow run main.nf \
  --csvFile input.csv \
  --run_busco_core \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro \
  --outdir results
```

### 5. Review Results

```bash
# Check BUSCO completeness
cat results/busco/*_short_summary.txt

# Check OMArk consistency
cat results/omark/*_detailed_summary.txt

# Review statistics
cat results/ensembl_stats/*.json
```

## Best Practices

!!! tip "Run All Workflows"
    For production annotations, always run BUSCO, OMArk, and Ensembl Stats together for comprehensive QC.

!!! tip "Validate Before Applying"
    Generate and review statistics before using `--apply_*` flags to load data into databases.

!!! tip "Use Specific Lineages"
    Choose the most specific BUSCO/OMArk lineage for your organism for best results.

!!! tip "Track Over Time"
    Keep statistics outputs in version control to monitor quality trends across releases.

!!! tip "Document Exceptions"
    Some species have genuine biological variations (gene losses, duplications) that affect scores—document these.

## Common Workflows by Role

### Annotation Team

```bash
# Complete annotation QC
nextflow run main.nf \
  --csvFile new_annotations.csv \
  --run_busco_core \
  --busco_mode both \
  --run_omark \
  --run_ensembl_stats \
  --host mysql-server.example.com \
  --user_r ensro \
  --enscode /software/ensembl/ENSCODE
```

### Assembly Team

```bash
# Assembly quality assessment
nextflow run main.nf \
  --csvFile assemblies.csv \
  --run_busco_ncbi \
  --outdir assembly_qc
```

### Database Administrator

```bash
# Generate and apply statistics
nextflow run main.nf \
  --csvFile release_dbs.csv \
  --run_ensembl_stats \
  --apply_ensembl_stats \
  --host mysql-server.example.com \
  --user ensadmin \
  --password ${DB_PASS} \
  --enscode /software/ensembl/ENSCODE \
  --team genebuild
```

### Comparative Genomics

```bash
# Multi-species comparison
nextflow run main.nf \
  --csvFile species_set.csv \
  --run_busco_core \
  --busco_mode protein \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro \
  --outdir comparative_qc
```
## Module Overview

The statistics pipeline consists of 13 modules organized into functional categories:

### Data Retrieval Modules
1. **[fetch-genome](modules/fetch-genome.md)** - Retrieves genome sequences from Ensembl core databases
2. **[fetch-proteins](modules/fetch-proteins.md)** - Extracts protein translations from Ensembl databases

### BUSCO Quality Assessment Modules
3. **[busco-dataset](modules/busco-dataset.md)** - Downloads appropriate BUSCO lineage datasets
4. **[busco-genome-lineage](modules/busco-genome-lineage.md)** - Runs BUSCO assessment on genome sequences
5. **[busco-protein-lineage](modules/busco-protein-lineage.md)** - Runs BUSCO assessment on protein translations
6. **[busco-core-metakeys](modules/busco-core-metakeys.md)** - Patches BUSCO metadata into core databases

### Orthology Analysis Modules
7. **[omamer-hog](modules/omamer-hog.md)** - Performs orthology inference using OMAmer
8. **[omark](modules/omark.md)** - Quality assessment of protein annotations using OMark

### Statistics Generation Modules
9. **[run-statistics](modules/run-statistics.md)** - Generates comprehensive annotation statistics
10. **[run-ensembl-meta](modules/run-ensembl-meta.md)** - Generates core database metadata SQL files

### Database Operations Modules
11. **[populate-db](modules/populate-db.md)** - Executes SQL files to populate databases
12. **[db-metadata](modules/db-metadata.md)** - Manages database metadata and versioning


## Pipeline Flow

The typical execution flow of the statistics pipeline:

```
1. Data Retrieval
   └─> FETCH_GENOME
   └─> FETCH_PROTEINS

2. Quality Assessment (Parallel)
   ├─> BUSCO_DATASET
   │   ├─> BUSCO_GENOME_LINEAGE
   │   └─> BUSCO_PROTEIN_LINEAGE
   │       └─> BUSCO_CORE_METAKEYS
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

6. Cleanup
   └─> CLEANING
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

## Support

### Documentation

- **Workflow Guides**: Detailed guides for each workflow
- **API Reference**: Parameter and configuration options
- **Troubleshooting**: Common issues and solutions
- **Examples**: Real-world use cases

### Getting Help

- **Issues**: [GitHub Issues](https://github.com/Ensembl/ensembl-genes/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Ensembl/ensembl-genes/discussions)
- **Contact**: Ensembl Genebuild team

## Related Documentation

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [BUSCO Documentation](https://busco.ezlab.org/busco_userguide.html)
- [OMArk Documentation](https://github.com/DessimozLab/OMArk)
- [Ensembl API Documentation](https://www.ensembl.org/info/docs/api/)

## Citation

If you use this pipeline, please cite:

```
Ensembl Genes Statistics Pipeline
https://github.com/Ensembl/ensembl-genes
```

And the relevant tools:
- **BUSCO**: Manni et al. (2021). DOI: 10.1093/molbev/msab199
- **OMArk**: Nevers et al. (2022). DOI: 10.1101/2022.11.25.517970
- **Ensembl**: Cunningham et al. (2022). DOI: 10.1093/nar/gkab1049


