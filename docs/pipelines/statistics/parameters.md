# Parameters Reference

Complete reference for all pipeline parameters.

## Required Parameters

### Core Input/Output

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--csvFile` | string | Path to input CSV metadata file | `data/genomes.csv` |
| `--outdir` | string | Output directory for results | `./results` |

!!! warning "Required"
    These parameters must be provided for every pipeline run.

## Workflow Selection

Enable specific analysis workflows:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_busco_core` | boolean | `false` | Run BUSCO using MySQL core database |
| `--run_busco_ncbi` | boolean | `false` | Run BUSCO using NCBI assembly accession (genome mode only) |
| `--run_omark` | boolean | `false` | Run OMArk proteome assessment |
| `--run_ensembl_stats` | boolean | `false` | Generate Ensembl statistics |
| `--run_ensembl_beta_metakeys` | boolean | `false` | Generate Ensembl beta metakeys |

### Examples

```bash
# Run BUSCO only
--run_busco_core

# Run multiple workflows
--run_busco_core --run_omark --run_ensembl_stats

# BUSCO from NCBI
--run_busco_ncbi
```

## Database Connection

Parameters for connecting to Ensembl core databases:

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--host` | string | Database host server | `mysql-server.example.com` |
| `--port` | integer | Database port | `3306` |
| `--user_r` | string | Read-only database user | `readonly_user` |
| `--user` | string | Database user with write permissions | `ensadmin` |
| `--password` | string | Database password | `secret123` |

!!! info "Database Access"
    - `--user_r` is sufficient for quality metrics generation
    - `--user` and `--password` are only required when applying statistics to the database

### Example

```bash
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_ensembl_stats \
  --host mysql-ens-sta-5.example.com \
  --port 4686 \
  --user_r ensro
```

## BUSCO Parameters

### Mode Selection

| Parameter | Type | Default | Options | Description |
|-----------|------|---------|---------|-------------|
| `--busco_mode` | string | `protein` | `protein`, `genome`, `both` | BUSCO analysis mode |

- **protein**: Analyze protein sequences from gene predictions
- **genome**: Analyze genome assembly directly
- **both**: Run both protein and genome modes

### Lineage Dataset

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--busco_dataset` | string | Default BUSCO lineage for all samples | `vertebrata_odb12` |

!!! tip "Per-Sample Lineage"
    You can also specify lineage per sample in the CSV file using the `busco_dataset` column.

### BUSCO Configuration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--busco_version` | string | `v6.0.0_cv1` | BUSCO container image version |
| `--download_path` | string | `/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/data/busco_data/data_odb12/` | Path to BUSCO lineage datasets |
| `--busco_datasets_file` | string | `../data/busco_lineage.json` | JSON file with available BUSCO lineages |
| `--dump_params` | string | `--canonical_only` | Parameters for sequence dumping |

### Available BUSCO Lineages

Common lineages (from OrthoDB v12):

| Lineage | Taxonomic Scope | Use For |
|---------|----------------|---------|
| `eukaryota_odb12` | All eukaryotes | Universal baseline |
| `metazoa_odb12` | Animals | All animal genomes |
| `vertebrata_odb12` | Vertebrates | Fish, amphibians, reptiles, birds, mammals |
| `mammalia_odb12` | Mammals | Human, mouse, cow, etc. |
| `primates_odb12` | Primates | Human, chimp, gorilla, etc. |
| `aves_odb12` | Birds | Chicken, zebra finch, etc. |
| `actinopterygii_odb12` | Ray-finned fish | Zebrafish, medaka, etc. |
| `insecta_odb12` | Insects | Fly, mosquito, bee, etc. |
| `diptera_odb12` | Flies | *Drosophila*, mosquitoes |
| `fungi_odb12` | Fungi | Yeast, *Aspergillus*, etc. |
| `viridiplantae_odb12` | Plants | All plant genomes |
| `embryophyta_odb12` | Land plants | Arabidopsis, rice, etc. |

For the complete list, check the BUSCO datasets file or visit [BUSCO website](https://busco.ezlab.org/).

### Example

```bash
# Use vertebrate lineage in protein mode
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --busco_mode protein \
  --busco_dataset vertebrata_odb12

# Run both protein and genome modes
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --busco_mode both
```

## OMArk Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--omamer_database` | string | `/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/data/omamer_db/LUCA_MinFamSize6_OR_MinFamComp05_A21_k6.h5` | Path to OMArk/OMAmer database |
| `--omark_singularity_path` | string | `/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/omark.sif` | Path to OMArk Singularity container |

### Example

```bash
nextflow run main.nf \
  --csvFile proteomes.csv \
  --run_omark \
  --omamer_database /data/omamer/LUCA.h5 \
  --host mysql-server.example.com \
  --user_r ensro
```

## Ensembl Statistics Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--enscode` | string | - | Path to Ensembl API/modules directory |
| `--bioperl` | string | `/bioperl-1.6.924` | Path to BioPerl installation |
| `--mysql_ensadmin` | string | `/hps/software/users/ensembl/ensw/mysql-cmds/ensembl/ensadmin` | Path to ensadmin script |
| `--meta_query_file` | string | `../bin/meta.sql` | SQL query file for metadata |
| `--project` | string | `ensembl` | Project name for metadata |
| `--team` | string | - | Team responsible (metakey) |

### Apply Statistics to Database

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--apply_ensembl_stats` | boolean | `false` | Insert statistics into the database |
| `--apply_ensembl_beta_metakeys` | boolean | `false` | Insert beta metakeys into the database |
| `--apply_busco_metakeys` | boolean | `false` | Create and load BUSCO metakeys JSON |

!!! warning "Database Write Access Required"
    When using `--apply_*` parameters, you must provide `--user` and `--password` with write permissions.

### Example

```bash
# Generate statistics only
nextflow run main.nf \
  --csvFile databases.csv \
  --run_ensembl_stats \
  --enscode /nfs/software/ensembl/ENSCODE \
  --host mysql-server.example.com \
  --user_r ensro

# Generate and apply to database
nextflow run main.nf \
  --csvFile databases.csv \
  --run_ensembl_stats \
  --apply_ensembl_stats \
  --enscode /nfs/software/ensembl/ENSCODE \
  --host mysql-server.example.com \
  --user ensadmin \
  --password secret123 \
  --team genebuild
```

## NCBI Download Parameters

For BUSCO NCBI mode:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--ncbiBaseUrl` | string | `https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/` | NCBI API base URL |

### Example

```bash
nextflow run main.nf \
  --csvFile ncbi_assemblies.csv \
  --run_busco_ncbi \
  --busco_dataset vertebrata_odb12
```

## Cache and Cleanup

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--cacheDir` | string | `/cache` | Directory for caching downloaded files |
| `--cleanCache` | boolean | `true` | Clean cache directory after pipeline completion |
| `--files_latency` | integer | `60` | File system latency in seconds |

### Example

```bash
# Keep cache for debugging
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_ncbi \
  --cleanCache false

# Use custom cache location
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_ncbi \
  --cacheDir /scratch/pipeline_cache
```

## Pipeline Info

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--tracedir` | string | `./results/pipeline_info` | Directory for pipeline execution reports |

The trace directory contains:

- `execution_trace.txt`: Task-level execution details
- `execution_timeline.html`: Visual timeline of task execution
- `execution_report.html`: Resource usage report
- `software_versions.yml`: Versions of all software used

## Resource Limits

Control maximum resources used by the pipeline:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--max_cpus` | integer | - | Maximum CPUs per process |
| `--max_memory` | memory | - | Maximum memory per process |
| `--max_time` | time | - | Maximum time per process |

### Example

```bash
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --max_cpus 16 \
  --max_memory 64.GB \
  --max_time 24.h
```

## Advanced Parameters

### Hidden Parameters

These are typically auto-configured but can be overridden:

| Parameter | Type | Description |
|-----------|------|-------------|
| `--dbname` | string | Database name (typically from CSV) |
| `--readme` | string | Path to README file |

## Complete Example

Full pipeline run with all common parameters:

```bash
nextflow run main.nf \
  --csvFile input_genomes.csv \
  --outdir /data/results/qc_pipeline \
  --run_busco_core \
  --run_omark \
  --run_ensembl_stats \
  --busco_mode both \
  --busco_dataset vertebrata_odb12 \
  --host mysql-ens-sta-5.ebi.ac.uk \
  --port 4686 \
  --user_r ensro \
  --enscode /nfs/software/ensembl/ENSCODE \
  --cacheDir /scratch/cache \
  --cleanCache true \
  --max_cpus 32 \
  --max_memory 128.GB \
  -profile singularity \
  -resume
```

## Parameter Files

For complex configurations, use a parameter file:

```yaml
# params.yml
csvFile: "genomes.csv"
outdir: "results"
run_busco_core: true
run_omark: true
busco_mode: "both"
host: "mysql-server.example.com"
port: 3306
user_r: "ensro"
enscode: "/software/ensembl/ENSCODE"
max_cpus: 32
max_memory: "128 GB"
```

Run with:

```bash
nextflow run main.nf -params-file params.yml -profile docker
```

## Next Steps

- [Input Format](input.md) - Learn how to structure your CSV input file
- [Output Documentation](output.md) - Understand the output files
- [Troubleshooting](troubleshooting.md) - Solve common parameter issues
