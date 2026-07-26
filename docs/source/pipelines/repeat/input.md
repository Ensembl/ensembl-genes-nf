# Input Format

The repeat pipeline uses a CSV file to specify the genomes to analyse and the metadata required to retrieve repeat libraries and annotate repetitive elements.

## CSV File Structure

### Basic Format

The CSV file must contain a header row followed by one row per genome.

```csv
column1,column2,column3
value1,value2,value3
value1,value2,value3
```

!!! important "CSV Requirements"
- **Header row required**: The first row must contain column names.
- **Comma-separated**: Use commas as delimiters.
- **One genome per row**: Each row represents a single assembly.
- **No trailing spaces**: Remove leading and trailing whitespace from all values.

---

## Required Columns

The repeat pipeline requires the following information for each genome.

| Column                  | Type    | Required    | Description                           | Example                                        |
| ----------------------- | ------- | ----------- | ------------------------------------- | ---------------------------------------------- |
| `gca`                   | string  | ✅           | NCBI assembly accession               | `GCA_000001405.29`                             |
| `taxon_id`              | integer | ✅           | NCBI taxonomy identifier              | `9606`                                         |
| `repeatmodeler_library` | string  | ⚠️ Optional | path to a RepeatModeler repeat library | `repeatmodeler.fa` |

If no RepeatModeler library is supplied, the pipeline can fall back to the configured repeat annotation strategy (for example Dfam or de novo repeat discovery, depending on the selected parameters).

---

## Example CSV

```csv
gca,taxon_id,repeatmodeler_library
GCA_000001405.29,9606,https://ftp.ensembl.org/pub/repeats/homo_sapiens.repeatmodeler.fa
GCA_000001635.9,10090,https://ftp.ensembl.org/pub/repeats/mus_musculus.repeatmodeler.fa
GCA_000002035.4,7955,
```

---

## Typical Command

```bash
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_repeatmasker \
  --run_red
```

---

# Column Details

## `gca`

NCBI GenBank genome assembly accession.

Format:

```text
GCA_000000000.0
```

Examples:

* `GCA_000001405.29`
* `GCA_000001635.9`
* `GCA_000002035.4`

Assemblies are downloaded automatically from NCBI when required.

---

## `taxon_id`

NCBI Taxonomy identifier.

Examples:

| Species                   | Taxonomy ID |
| ------------------------- | ----------: |
| *Homo sapiens*            |        9606 |
| *Mus musculus*            |       10090 |
| *Danio rerio*             |        7955 |
| *Drosophila melanogaster* |        7227 |

Taxonomy identifiers are used to select appropriate repeat resources when available.

---

## `repeatmodeler_library`

Optional path pointing to a RepeatModeler repeat library.

Example:

```text
path/homo_sapiens.repeatmodeler.fa
```

If omitted, the pipeline uses the configured fallback strategy.

---

# Complete Example

```csv
gca,taxon_id,repeatmodeler_library
GCA_000001405.29,9606,path/homo_sapiens.repeatmodeler.fa
GCA_000001635.9,10090,path/mus_musculus.repeatmodeler.fa
GCA_000002035.4,7955,
```

Command:

```bash
nextflow run main.nf \
    --csvFile genomes.csv \
    --run_repeatmasker \
    --run_red \
    --outdir results
```

---

# Validation

The input CSV is validated before the pipeline starts.

Typical validation errors include:

### ❌ Invalid GCA accession

```text
Error: GCA must follow the format GCA_000000000.0
```

### ❌ Missing taxonomy identifier

```text
Error: taxon_id is required
```

### ❌ Invalid RepeatModeler path

```text
Error: RepeatMasker library not specified
```

The `CHECK_AND_DOWNLOAD_RMLIBRARY` process validates remote URLs before the workflow continues.

---

# Tips and Best Practices

!!! tip "Use public RepeatModeler libraries"

Whenever possible, provide curated RepeatModeler libraries generated for the target assembly.

!!! tip "Verify URLs"

Ensure URLs are accessible before launching the pipeline.

```bash
wget --spider https://your.server/repeatmodeler.fa
```

!!! tip "Use stable GCA accessions"

Always use the final assembly accession rather than intermediate versions.

!!! tip "Run RED together with RepeatMasker"

Using both methods provides complementary repeat annotations.

---

# Schema Validation

The CSV file is validated against the repeat pipeline input schema.

```bash
cat pipelines/repeat/assets/schema_input.json
```

You can also validate the inputs by running:

```bash
nextflow run main.nf --csvFile genomes.csv --help
```

---

# Next Steps

* [Parameters Reference](parameters.md)
* [Output Documentation](output.md)
* [RepeatMasker Workflow](workflows/index.md)
* [Module Documentation](modules/index.md)
