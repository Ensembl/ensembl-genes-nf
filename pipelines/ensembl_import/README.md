# Ensembl Import Pipeline

This pipeline imports RefSeq assemblies into Ensembl core databases. For each
assembly it fetches the RefSeq annotation, creates a core database, adds sample
gene metadata, loads taxonomy information from NCBI, and publishes an input
CSV for the statistics pipeline.

## Requirements

- Nextflow `26.04.0` or a compatible 26.x release
- Singularity for normal execution
- A local checkout of `ensembl-genes` containing:
  - `src/python/ensembl/genes/ensembl_loading/gff_cli.py`
  - `src/python/ensembl/genes/metadata/core_metadata.py`
- Network access to RefSeq and NCBI for non-stub runs
- MySQL access to the target server

The pipeline uses the `nf-schema` plugin to validate command-line parameters
and the sample sheet before work starts.

## Sample Sheet

Pass the sample sheet with `--input_csv`. It must be a CSV with a header and
one or more rows:

```csv
gcf,species
GCF_000001405.40,Homo sapiens
GCA_900184115.1,Example species
```

The `gcf` value must be a RefSeq accession in the form `GCA_#########.#` or
`GCF_#########.##`. The `species` value must not be empty.

The row schema is defined in
[`assets/schema_input.json`](assets/schema_input.json).

## Database Settings

Pass a JSON file with `--server_settings`:

```json
{
  "db_host": "mysql.example.org",
  "db_port": 3306,
  "db_user": "writeuser",
  "db_password": "secret",
  "db_read_user": "readuser"
}
```

`db_user` is used for database creation and metadata loading. `db_read_user`
is used by the taxonomy step to read the taxonomy ID from the core database.

## Running

The pipeline entrypoint is `main.nf`. From the repository root:

```bash
nextflow run pipelines/ensembl_import/main.nf \
  --input_csv /path/to/assemblies.csv \
  --server_settings /path/to/server_settings.json \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --outdir /path/to/results
```

The default executor is Slurm. For a local stub check:

```bash
nextflow run pipelines/ensembl_import/main.nf \
  -profile local \
  -stub-run \
  --input_csv /path/to/assemblies.csv \
  --server_settings /path/to/server_settings.json \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --outdir /tmp/ensembl_import_results
```

The stub run validates the workflow graph and channel wiring without running
the database, RefSeq, or NCBI commands.

## Workflow

The workflow is structured as:

```text
main.nf
  workflows/import_refseq.nf
    subworkflows/import_refseq_to_core.nf
      FETCH_REFSEQ
      LOAD_REFSEQ
    subworkflows/get_metadata_core.nf
      GET_SAMPLE_GENE
      ADD_STATIC_METAKEYS
      LOAD_TAXONOMY
    PREPARE_STATS_INPUT
```

Each assembly is processed independently.

1. `FETCH_REFSEQ` downloads and prepares the RefSeq GFF3, FASTA, and assembly
   report files.
2. `LOAD_REFSEQ` creates the Ensembl core database with `gff_cli.py`.
3. `GET_SAMPLE_GENE` selects and inserts sample gene metadata.
4. `ADD_STATIC_METAKEYS` inserts configured static metadata keys.
5. `LOAD_TAXONOMY` retrieves taxonomy data from NCBI and inserts it into the
  core database.
6. `PREPARE_STATS_INPUT` writes the CSV consumed by the statistics pipeline.

## Outputs

Published files are written below `--outdir`:

```text
<outdir>/refseq/<assembly_id>/
<outdir>/metadata/<assembly_id>/
```

The workflow emits:

- RefSeq input files and process completion markers
- Metadata and taxonomy completion markers
- `statistics_input.csv`, containing the core name, species ID, and the
  published genome FASTA path
- `versions.yml` files from each process

The generated CSV uses `UNKNOWN` for `taxon_id` and `gca`, allowing the
statistics pipeline to query those values from the loaded core. The
`genome_file` column points to the published RefSeq FASTA, so the statistics
pipeline does not download the genome again. No protein file is supplied;
protein sequences can be dumped from the core by the statistics pipeline.

## Configuration

Pipeline-specific configuration is in
[`nextflow.config`](nextflow.config). The default profile uses Slurm. Use the
`local` profile for development and stub runs.

The parameter schema is in
[`nextflow_schema.json`](nextflow_schema.json), and the sample-sheet schema is
in [`assets/schema_input.json`](assets/schema_input.json).

The registry update modules are not part of the runnable import workflow yet.
