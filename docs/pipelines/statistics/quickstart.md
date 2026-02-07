# Quick Start Guide

Get started with the Ensembl Genes Statistics pipeline in minutes!

## Prerequisites

Before you begin, ensure you have:

- [x] **Nextflow** (version 23.04 or later)
- [x] **Java** (version 11 or later)
- [x] **Container engine** (Docker, Singularity, or Podman)
- [x] Access to genome data or Ensembl core databases

## Installation

### 1. Install Nextflow

If you haven't already installed Nextflow:

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

Verify installation:

```bash
nextflow -version
```

### 2. Clone the Repository

```bash
git clone https://github.com/Ensembl/ensembl-genes-nf.git
cd ensembl-genes-nf/pipelines/statistics
```

## Basic Usage

### Example 1: BUSCO Analysis with NCBI Assembly

This is the simplest way to get started—analyze a genome directly from NCBI:

```bash
# Create input CSV
cat > my_genome.csv << EOF
gca,taxon_id,busco_dataset
GCA_000001405.29,9606,primates_odb12
EOF

# Run the pipeline
nextflow run main.nf \
  --csvFile my_genome.csv \
  --run_busco_ncbi \
  --outdir results \
  -profile docker
```

!!! success "What happens next?"
    The pipeline will:
    
    1. Download the genome assembly from NCBI
    2. Download the appropriate BUSCO lineage dataset
    3. Run BUSCO in genome mode
    4. Generate quality reports in `results/busco/`

### Example 2: BUSCO Analysis from Core Database

If you have an Ensembl core database:

```bash
# Create input CSV
cat > my_genomes.csv << EOF
dbname,species_id,busco_dataset,taxon_id
homo_sapiens_core_110_38,1,primates_odb12,9606
mus_musculus_core_110_39,1,glires_odb12,10090
EOF

# Run the pipeline
nextflow run main.nf \
  --csvFile my_genomes.csv \
  --run_busco_core \
  --busco_mode both \
  --host mysql-server.example.com \
  --port 3306 \
  --user_r readonly_user \
  --outdir results \
  -profile singularity
```

!!! info "Database Access"
    You need read-only access to the MySQL server hosting the core databases.

### Example 3: OMArk Proteome Assessment

```bash
cat > proteomes.csv << EOF
dbname,species_id,taxon_id
danio_rerio_core_110_11,1,7955
xenopus_tropicalis_core_110_10,1,8364
EOF

nextflow run main.nf \
  --csvFile proteomes.csv \
  --run_omark \
  --host mysql-server.example.com \
  --port 3306 \
  --user_r readonly_user \
  --outdir results \
  -profile docker
```

### Example 4: Ensembl Statistics Generation

```bash
cat > databases.csv << EOF
dbname,species_id
caenorhabditis_elegans_core_110_280,1
drosophila_melanogaster_core_110_9,1
EOF

nextflow run main.nf \
  --csvFile databases.csv \
  --run_ensembl_stats \
  --host mysql-server.example.com \
  --port 3306 \
  --user_r readonly_user \
  --enscode /path/to/ensembl/modules \
  --outdir results \
  -profile singularity
```

### Example 5: Complete Quality Control Pipeline

Run all quality metrics together:

```bash
cat > complete_qc.csv << EOF
dbname,species_id,busco_dataset,taxon_id,protein_file,genome_file
arabidopsis_thaliana_core_110_11,1,brassicales_odb12,3702,/data/proteins.fa,/data/genome.fa
EOF

nextflow run main.nf \
  --csvFile complete_qc.csv \
  --run_busco_core \
  --run_omark \
  --run_ensembl_stats \
  --busco_mode both \
  --host mysql-server.example.com \
  --port 3306 \
  --user_r readonly_user \
  --enscode /path/to/ensembl/modules \
  --outdir results \
  -profile docker
```

## Configuration Profiles

The pipeline supports several execution profiles:

| Profile | Description | Command |
|---------|-------------|---------|
| `docker` | Use Docker containers | `-profile docker` |
| `singularity` | Use Singularity containers | `-profile singularity` |
| `conda` | Use Conda environments | `-profile conda` |
| `test` | Run with test data | `-profile test` |

You can combine profiles:

```bash
nextflow run main.nf -profile docker,test
```

## Understanding Results

After the pipeline completes, check your results:

```bash
# View directory structure
tree results/

# BUSCO summary
cat results/busco/*_busco_short_summary.txt

# OMArk summary
cat results/omark/*_omark_proteins_detailed_summary.txt

# Statistics JSON
cat results/ensembl_stats/*_statistics.json
```

## Common Options

### Specify Output Directory

```bash
--outdir /path/to/output
```

### Resume Failed Run

```bash
nextflow run main.nf -resume --csvFile my_genomes.csv --run_busco_core
```

!!! tip "Always use -resume"
    The `-resume` flag allows Nextflow to skip completed tasks, saving time and compute resources.

### Run in Background

```bash
nextflow run main.nf --csvFile genomes.csv --run_busco_core -bg > pipeline.log 2>&1
```

### Custom BUSCO Lineage

```bash
--busco_dataset vertebrata_odb12
```

Available lineages are automatically downloaded. See [BUSCO documentation](workflows/busco.md) for the full list.

### Clean Cache After Completion

```bash
--cleanCache true
```

This removes temporary files after successful completion.

## Resource Configuration

For large genomes or multiple samples, adjust resources:

```bash
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  -profile docker \
  --max_cpus 32 \
  --max_memory 128.GB \
  --max_time 72.h
```

## Monitoring Progress

### View Real-time Log

```bash
tail -f .nextflow.log
```

### Check Execution Report

After completion, Nextflow generates reports:

```bash
# Timeline visualization
firefox results/pipeline_info/execution_timeline.html

# Resource usage report
firefox results/pipeline_info/execution_report.html
```

## Troubleshooting Quick Tips

!!! warning "Pipeline fails immediately?"
    - Check input CSV format: [Input Format Guide](input.md)
    - Verify database connectivity: `mysql -h <host> -P <port> -u <user> -p`
    - Ensure container engine is running: `docker ps` or `singularity --version`

!!! warning "Out of memory errors?"
    - Increase memory limits: `--max_memory 256.GB`
    - Reduce parallel processes: `--max_cpus 8`

!!! warning "BUSCO lineage not found?"
    - Check available lineages: `--busco_datasets_file data/busco_lineage.json`
    - Verify download path: `--download_path /path/to/busco/data`

For more detailed troubleshooting, see the [Troubleshooting Guide](troubleshooting.md).

## Next Steps

- [Detailed Parameters](parameters.md) - Explore all available parameters
- [Input Format](input.md) - Learn about input file requirements
- [Output Documentation](output.md) - Understand the output files
- [BUSCO Workflow](workflows/busco.md) - Deep dive into BUSCO analysis
- [OMArk Workflow](workflows/omark.md) - Learn about proteome assessment
- [Ensembl Stats](workflows/ensembl-stats.md) - Statistics generation details

## Getting Help

If you encounter issues:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Review [GitHub Issues](https://github.com/Ensembl/ensembl-genes-nf/issues)
3. Contact [Ensembl Helpdesk](http://www.ensembl.org/Help/Contact)
