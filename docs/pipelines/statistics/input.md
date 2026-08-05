# Input Format

The pipeline uses a CSV file to specify samples and their metadata. This guide explains the required format and provides examples for different use cases.

## CSV File Structure

### Basic Format

The CSV file must contain a header row with column names, followed by one row per genome to analyze.

```csv
column1,column2,column3
value1,value2,value3
value1,value2,value3
```

!!! important "CSV Requirements"
    - **Header row required**: First row must contain column names
    - **Comma-separated**: Use commas as delimiters
    - **No trailing spaces**: Column values should not have leading/trailing whitespace
    - **One genome per row**: Each row represents a single sample

## Column Specifications by Workflow

Different workflows require different columns. Choose the columns based on which workflows you're running.

### For BUSCO (Core Database Mode)

Required when using `--run_busco_core`:

| Column | Type | Required | Description | Example |
|--------|------|----------|-------------|---------|
| `dbname` | string | ✅ | Ensembl core database name | `homo_sapiens_core_110_38` |
| `species_id` | integer | ✅ | Species ID in the core database | `1` |
| `busco_dataset` | string | ⚠️ * | Specific BUSCO lineage dataset (if not will select the closest one in terms of taxonomy ) | `primates_odb12` |
| `protein_file` | string | ⚠️ | Path to protein FASTA file (if not using DB) | `/data/proteins.fa` |
| `genome_file` | string | ⚠️ | Path to genome FASTA file (if not using DB) | `/data/genome.fa` |

\* Can be specified globally with `--busco_dataset` parameter instead

**Example CSV:**
```csv
dbname,species_id,busco_dataset
homo_sapiens_core_110_38,1,primates_odb12
mus_musculus_core_110_39,1,glires_odb12
danio_rerio_core_110_11,1,actinopterygii_odb12
```

### For BUSCO (NCBI Mode)

Required when using `--run_busco_ncbi`:

| Column | Type | Required | Description | Example |
|--------|------|----------|-------------|---------|
| `gca` | string | ✅ | NCBI assembly accession | `GCA_000001405.29` |
| `taxon_id` | integer | ✅ | NCBI taxonomy ID | `9606` |
| `busco_dataset` | string | ⚠️ * | Specific BUSCO lineage dataset (if not will select the closest one in terms of taxonomy ) | `primates_odb12` |

\* Can be specified globally with `--busco_dataset` parameter instead

**Example CSV:**
```csv
gca,taxon_id,busco_dataset
GCA_000001405.29,9606,primates_odb12
GCA_000001635.9,10090,glires_odb12
GCA_000002035.4,7955,actinopterygii_odb12
```

!!! info "NCBI Mode"
    NCBI mode only runs BUSCO in genome mode (not protein mode). The genome assembly is automatically downloaded from NCBI.

### For OMArk

Required when using `--run_omark`:

| Column | Type | Required | Description | Example |
|--------|------|----------|-------------|---------|
| `dbname` | string | ✅ | Ensembl core database name | `xenopus_tropicalis_core_110_10` |
| `species_id` | integer | ✅ | Species ID in the core database | `1` |
| `protein_file` | string | ⚠️ | Path to protein FASTA file (if not using DB) | `/data/proteins.fa` |

**Example CSV:**
```csv
dbname,species_id
xenopus_tropicalis_core_110_10,1
danio_rerio_core_110_11,1
```

### For Ensembl Statistics

Required when using `--run_ensembl_stats` or `--run_ensembl_beta_metakeys`:

| Column | Type | Required | Description | Example |
|--------|------|----------|-------------|---------|
| `dbname` | string | ✅ | Ensembl core database name | `arabidopsis_thaliana_core_110_11` |
| `species_id` | integer | ✅ | Species ID in the core database | `1` |

**Example CSV:**
```csv
dbname,species_id
arabidopsis_thaliana_core_110_11,1
caenorhabditis_elegans_core_110_280,1
drosophila_melanogaster_core_110_9,1
```

## Complete Examples

### Example 1: BUSCO Quality Control (Core DB)

Full quality assessment with both protein and genome BUSCO modes:

```csv
dbname,species_id,busco_dataset
homo_sapiens_core_110_38,1,primates_odb12
pan_troglodytes_core_110_40,1,primates_odb12
gorilla_gorilla_core_110_6,1,primates_odb12
pongo_abelii_core_110_5,1,primates_odb12
```

Command:
```bash
nextflow run main.nf \
  --csvFile primates.csv \
  --run_busco_core \
  --busco_mode both \
  --host mysql-server.example.com \
  --user_r ensro
```

### Example 2: BUSCO from NCBI Assemblies

Quick assessment of public genome assemblies:

```csv
gca,taxon_id,busco_dataset
GCA_000001405.29,9606,primates_odb12
GCA_011100555.1,9598,primates_odb12
GCA_008122165.1,9593,primates_odb12
```

Command:
```bash
nextflow run main.nf \
  --csvFile ncbi_primates.csv \
  --run_busco_ncbi
```

### Example 3: Complete Multi-Workflow Analysis

Run all quality metrics together:

```csv
dbname,species_id,busco_dataset
drosophila_melanogaster_core_110_9,1,diptera_odb12
anopheles_gambiae_core_110_56,1,diptera_odb12
aedes_aegypti_core_110_6,1,diptera_odb12
```

Command:
```bash
nextflow run main.nf \
  --csvFile insects.csv \
  --run_busco_core \
  --run_omark \
  --run_ensembl_stats \
  --busco_mode both \
  --host mysql-server.example.com \
  --user_r ensro \
  --legacy_enscode /software/ensembl/ENSCODE
```

### Example 4: Mixed Lineages

Different BUSCO lineages for different species:

```csv
dbname,species_id,busco_dataset
homo_sapiens_core_110_38,1,primates_odb12
danio_rerio_core_110_11,1,actinopterygii_odb12
drosophila_melanogaster_core_110_9,1,diptera_odb12
arabidopsis_thaliana_core_110_11,1,embryophyta_odb12
saccharomyces_cerevisiae_core_110_4,1,saccharomycetes_odb12
```

Command:
```bash
nextflow run main.nf \
  --csvFile diverse_species.csv \
  --run_busco_core \
  --host mysql-server.example.com \
  --user_r ensro
```

### Example 5: Using External FASTA Files

When sequences are not in a database:

```csv
dbname,species_id,busco_dataset,protein_file,genome_file
custom_genome_v1,1,vertebrata_odb12,/data/project1/proteins.fa,/data/project1/genome.fa
custom_genome_v2,1,vertebrata_odb12,/data/project2/proteins.fa,/data/project2/genome.fa
```

Command:
```bash
nextflow run main.nf \
  --csvFile custom_genomes.csv \
  --run_busco_core \
  --busco_mode both
```

## Column Details

### `dbname`

Format: `{species}_{type}_{release}_{assembly_version}`

- **species**: Species name (lowercase, underscores)
- **type**: Usually `core`
- **release**: Ensembl release number
- **assembly_version**: Assembly version number

Examples:
- `homo_sapiens_core_110_38` (Human, release 110, GRCh38)
- `mus_musculus_core_110_39` (Mouse, release 110, GRCm39)

### `species_id`

The species ID within the core database. In most cases, this is `1`.

!!! info "Multi-species databases"
    For bacteria collections or multi-species databases, this might be different. Check the `meta` table in your core database.

### `busco_dataset`

The BUSCO lineage dataset name. Must match available datasets in OrthoDB v12.

Format: `{lineage}_odb12`

Common lineages:
- `primates_odb12`, `mammalia_odb12`, `vertebrata_odb12`
- `actinopterygii_odb12`, `aves_odb12`
- `diptera_odb12`, `insecta_odb12`
- `embryophyta_odb12`, `viridiplantae_odb12`
- `fungi_odb12`, `saccharomycetes_odb12`

See [BUSCO documentation](workflows/busco.md) for complete list.

### `taxon_id`

The NCBI Taxonomy ID for the species.

Find taxonomy IDs at: [https://www.ncbi.nlm.nih.gov/taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

Examples:
- `9606` = *Homo sapiens* (Human)
- `10090` = *Mus musculus* (Mouse)
- `7955` = *Danio rerio* (Zebrafish)
- `7227` = *Drosophila melanogaster* (Fruit fly)
- `3702` = *Arabidopsis thaliana* (Thale cress)

### `gca`

NCBI GenBank assembly accession.

Format: `GCA_{9digits}.{version}`

Examples:
- `GCA_000001405.29` (Human GRCh38.p13)
- `GCA_000001635.9` (Mouse GRCm39)

Find assemblies at: [https://www.ncbi.nlm.nih.gov/assembly](https://www.ncbi.nlm.nih.gov/assembly)

### `protein_file` and `genome_file`

Full paths to FASTA files:

- **protein_file**: Protein sequences (canonical transcripts preferred)
- **genome_file**: Genome assembly (unmasked or soft-masked)

Requirements:
- Files must exist and be readable
- FASTA format (.fa, .fasta, .fna, .faa)
- Can be gzip-compressed (.gz)

## Validation

The pipeline validates your CSV file against a JSON schema. Common errors:

### ❌ Invalid GCA format
```
Error: GCA must be in the format GCA_000000000.0
```
**Fix**: Ensure GCA follows the pattern `GCA_{9digits}.{version}`

### ❌ Missing required column
```
Error: Column 'dbname' is required but not found
```
**Fix**: Add the missing column to your CSV header

### ❌ Invalid taxon_id
```
Error: Taxon ID must be a number
```
**Fix**: Use numeric NCBI taxonomy IDs only

### ❌ Trailing spaces
```
Error: Database names cannot contain or have trailing spaces
```
**Fix**: Remove spaces: `homo_sapiens_core_110_38` not `homo_sapiens_core_110_38 `

## Tips and Best Practices

!!! tip "Start small"
    Test with 1-2 samples first, then scale up to your full dataset.

!!! tip "Use absolute paths"
    For file paths (protein_file, genome_file), use absolute paths to avoid issues with working directories.

!!! tip "Check database connectivity first"
    Before running the pipeline, verify you can connect to the database:
    ```bash
    mysql -h mysql-server.example.com -P 3306 -u ensro -e "SHOW DATABASES LIKE 'homo_sapiens_core_%';"
    ```

!!! tip "Validate GCA accessions"
    Verify NCBI accessions exist before running:
    ```bash
    curl "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_000001405.29"
    ```

!!! tip "Group by lineage"
    For better resource efficiency, group species with the same BUSCO lineage together.

## Schema Validation

The input CSV is validated against a JSON schema. You can view the schema:

```bash
cat pipelines/statistics/assets/schema_input.json
```

Or validate manually:

```bash
# Using nf-schema plugin
nextflow run main.nf --csvFile genomes.csv --help
```

## Next Steps

- [Parameters Reference](parameters.md) - Configure pipeline behavior
- [Output Documentation](output.md) - Understand results
