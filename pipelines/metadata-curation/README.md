# Multi-Source Metadata Curation Pipeline

A Nextflow pipeline for extracting, standardizing, and curating metadata from multiple genomics data sources including SRA/ENA, GEO, and PMC publications.

## Quick Start

### Prerequisites

- Nextflow ≥ 23.04.0
- Mamba/Conda environment with Python 3.11 and requests library

### Setup Environment

```bash
# Create conda environment
mamba create -n metadata-curation python=3.11 requests -y

# Or use existing environment
mamba activate metadata-curation
```

### Run the Pipeline

```bash
# Basic usage
nextflow run main.nf --query "ribosome profiling" --outdir results

# With custom parameters
nextflow run main.nf \
    --query "RNA-seq AND cancer" \
    --max_results 50 \
    --outdir results \
    --output_format csv

# Using parameter file
nextflow run main.nf -params-file params.yaml
```

## Pipeline Overview

The pipeline implements a multi-stage metadata curation workflow:

1. **Data Source Extraction**
   - SRA/ENA metadata via NCBI E-utilities
   - GEO series and sample metadata
   - PMC publication content (optional)

2. **Ontology Standardization**
   - Term mapping to biomedical ontologies (ChEBI, CL, EFO, etc.)
   - Manual synonym handling for common terms
   - Confidence scoring for mappings

3. **Context-Aware Processing**
   - Study-level experimental design extraction
   - Sample-level metadata processing with sequencing-type schemas
   - Cross-modal biogroup detection

4. **Quality Validation & Output**
   - Metadata completeness validation
   - Biogroup quality assessment
   - Tabular output generation (CSV/TSV)

## Output Files

The pipeline generates:

- `sample_table.csv` - Sample-level metadata table
- `biogroups.json` - Detected biogroups with cross-modal analysis
- `validation_report.json` - Quality assessment results
- `ontology_cache/` - Cached ontology terms for reuse

## Key Features

### Multi-Source Integration
- Combines SRA experimental metadata with GEO annotations
- Links publication context from PMC full-text articles
- Hierarchical data source prioritization

### Sequencing-Type Awareness
- RNA-seq: Cell type, treatment, timepoint, replicate
- Ribo-seq: Translation inhibitors, fractionation methods
- CAGE: Cap selection methods, TSS enrichment
- Long-read: Platform specifics, library preparation

### Cross-Modal Biogroup Detection
- Groups samples by biological conditions across sequencing types
- Enables pairing of RNA-seq with Ribo-seq from same experiments
- Semantic equivalence mapping (e.g., "untreated" = "control")

### Quality Assessment
- Metadata completeness scoring per sequencing type
- Ontology mapping success rates
- Biogroup detection quality metrics

## Configuration

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `query` | - | Search query for studies (required) |
| `max_results` | 1000 | Maximum studies per source |
| `include_pmc` | false | Include PMC publication extraction |
| `email` | null | Email for NCBI API requests |
| `output_format` | csv | Output format (csv/tsv/json) |
| `output_level` | sample | Output granularity (sample/biogroup/study) |
| `ontology_cache_dir` | ontology_cache | Cache directory for ontology terms |

### Profiles

- `test` - Minimal test dataset (3 studies)
- `mamba` - Use mamba for dependency management
- `docker` - Use containerized execution

## Example Queries

```bash
# Ribosome profiling studies
nextflow run test_simple.nf --query "ribosome profiling"

# Cancer RNA-seq studies  
nextflow run test_simple.nf --query "RNA-seq AND cancer"

# Cell type specific studies
nextflow run test_simple.nf --query "HEK293 AND transcriptome"

# Time course experiments
nextflow run test_simple.nf --query "time course AND RNA-seq"
```

## Output Schema

### Sample Table Columns

| Column | Description |
|--------|-------------|
| sample_id | Run/sample accession |
| experiment_id | Experiment accession |
| study_id | BioProject/GSE ID |
| sequencing_type | Detected sequencing technology |
| cell_type | Standardized cell type |
| treatment | Treatment/condition |
| timepoint | Time point (if applicable) |
| organism | Organism name |
| platform | Sequencing platform |
| total_reads | Number of reads |

### Biogroup Schema

Biogroups represent sets of samples with identical biological conditions that can be compared across sequencing modalities.

## Implementation Notes

- **Rate Limiting**: Respects NCBI API rate limits (3 requests/second)
- **Caching**: Ontology terms cached locally for performance
- **Error Handling**: Robust error handling with detailed logging
- **Scalability**: Modular design supports large-scale processing
- **Transparency**: Preserves original metadata with mapping provenance

## Testing

The pipeline has been tested with real data:
- SRA extraction: 3 studies, 46 experiments, 47 runs
- Sample processing: Cell type and treatment extraction
- Table generation: 47 sample rows with standardized metadata

Example output shows successful extraction of CFTR overexpression studies with proper cell type detection (HEK293 cells) and treatment parsing.