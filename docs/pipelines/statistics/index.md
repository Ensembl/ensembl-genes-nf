# Statistics Pipeline

The **Statistics Pipeline** generates comprehensive quality metrics and statistics for gene annotations and assemblies. This pipeline is essential for validating annotations, assessing completeness, and providing metadata for Ensembl databases.

## Quick Links

- **[Quick Start Guide](quickstart.md)** - Get up and running quickly
- **[Input Specification](input.md)** - Prepare your input data
- **[Output Reference](output.md)** - Understand the results
- **[Configuration Guide](configuration.md)** - Customize the pipeline
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
