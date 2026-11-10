# Ensembl Import Pipeline

This pipeline imports RefSeq assemblies into Ensembl core databases. For each
assembly it fetches the RefSeq annotation, creates a core database, optionally
loads RefSeq RepeatMasker features, adds sample gene and static metadata, loads
taxonomy information from NCBI, updates the registry, and publishes input CSVs
for downstream processing.

## Requirements

- Nextflow `26.04.0` or a compatible 26.x release
- Singularity for normal execution
- A local installation of `ensembl-genes` (feature/refseq_loading_anna)
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
```

The `gcf` value must be a RefSeq accession in the form `GCF_#########.##`. The
`species` value must be a binomial name consisting of
two words.

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
    "db_read_user": "readuser",
    "registry_host": "mysql.registry.example.org",
    "registry_port": 4527,
    "registry_user": "registryuser",
    "registry_password": "registry-secret",
    "registry_database": "refseq_assembly"
}
```

`db_user` is used for database creation and metadata updates. `db_read_user`
is used by the taxonomy step to read the taxonomy ID from the core database.
The `registry_*` fields are used by the registry tracker. Registry connection
details are required; the pipeline has no hardcoded registry fallback values.

## Running

The pipeline entrypoint is `main.nf`. From the repository root:

```bash
nextflow run pipelines/ensembl_import/main.nf \
  --input_csv /path/to/assemblies.csv \
  --server_settings /path/to/server_settings.json \
  --genebuilder automated \
  --skip_start_check false \
  --outdir /path/to/results
```

The default executor is Slurm. For a local stub check:

```bash
nextflow run pipelines/ensembl_import/main.nf \
  -profile local \
  -stub-run \
  --input_csv /path/to/assemblies.csv \
  --server_settings /path/to/server_settings.json \
  --outdir /tmp/ensembl_import_results
```

The stub run validates the workflow graph and channel wiring without running
the database, RefSeq, or NCBI commands.

## Workflow

The workflow is structured as:

```text
main.nf
  workflows/import_refseq.nf
    REGISTRY_TRACKER_START
    subworkflows/import_refseq_to_core.nf
      FETCH_REFSEQ
      LOAD_REFSEQ
    subworkflows/import_repeatmasker.nf
      DOWNLOAD_REPEATMASKER
      CONVERT_REPEATMASKER
      LOAD_REPEATMASKER
    subworkflows/get_metadata_core.nf
      GET_SAMPLE_GENE
      ADD_STATIC_METAKEYS
      LOAD_TAXONOMY
    GET_CORE_METADATA
    PREPARE_CORE_RENAME
    REGISTRY_TRACKER_COMPLETE
    PREPARE_STATS_INPUT
    COMBINE_CORE_RENAME
    COMBINE_STATS_INPUT
```

Each assembly is processed independently.

1. `REGISTRY_TRACKER_START` checks the registry and records the import as
   `in_progress`. The pipeline stops if the current GCF already has an
   import-status record.
2. `FETCH_REFSEQ` downloads and prepares the RefSeq GFF3, FASTA, and assembly
   report files.
3. `LOAD_REFSEQ` creates the Ensembl core database with `gff_cli.py`.
4. `DOWNLOAD_REPEATMASKER` attempts to download the RefSeq RepeatMasker output
   with `gff-loader refseq download-repeatmasker`. A 404 is recorded as a
   missing optional file and does not fail the import.
5. `CONVERT_REPEATMASKER` converts the downloaded output to GTF using the
   assembly report when the file is available.
6. `LOAD_REPEATMASKER` loads the GTF with analysis name `repeatmasker` when
   conversion produced one.
7. `GET_SAMPLE_GENE` selects and inserts sample gene metadata.
8. `ADD_STATIC_METAKEYS` inserts configured static metadata keys.
9. `LOAD_TAXONOMY` retrieves taxonomy data from NCBI and inserts it into the
   core database. Requests are serialized and transient NCBI rate-limit or
   network errors are retried with backoff.
10. `GET_CORE_METADATA` reads `genebuild.version` and
    `genebuild.last_geneset_update` from each imported core.
11. `PREPARE_CORE_RENAME` creates a per-assembly mapping from the current core
    name to the next core name, using the update month as a `YYYYMM` suffix.
12. `REGISTRY_TRACKER_COMPLETE` writes the core genebuild version and update
    date to the registry and marks the import as `completed`.
13. `PREPARE_STATS_INPUT` writes the per-assembly CSV consumed by the
    statistics pipeline.
14. `COMBINE_CORE_RENAME` and `COMBINE_STATS_INPUT` combine the per-assembly
    CSV files after all assemblies have reached the final output stage.

## Outputs

Published files are written below `--outdir`:

```text
<outdir>/refseq/<assembly_id>/
<outdir>/refseq/<assembly_id>/repeatmasker/
<outdir>/metadata/<assembly_id>/
<outdir>/core_rename/core_rename.csv
<outdir>/statistics_input/statistics_input.csv
```

The workflow emits:

- RefSeq input files and process completion markers
- Metadata and taxonomy completion markers
- `statistics_input.csv`, containing the core name, species ID, and the
  published genome FASTA path
- `core_rename.csv`, containing `core_name`, `server`, and `new_core_name`
- `versions.yml` files from each process

The `core_rename.csv` output uses the core’s `genebuild.last_geneset_update`
value. For example, a core ending in `_1` and an update month of `2026-01`
produces a new name ending in `_202601`.

The generated CSV contains only `dbname`, `species_id`, and `genome_file`.
The `genome_file` column points to the published RefSeq FASTA, so the
statistics pipeline does not download the genome again. No protein file is
supplied; protein sequences can be dumped from the core by the statistics
pipeline.

## Configuration

Pipeline-specific configuration is in
[`nextflow.config`](nextflow.config). The default profile uses Slurm. Use the
`local` profile for development and stub runs.

The parameter schema is in
[`nextflow_schema.json`](nextflow_schema.json), and the sample-sheet schema is
in [`assets/schema_input.json`](assets/schema_input.json).

The registry tracker is integrated into the workflow. It can also be called
directly for inspection or recovery:

```bash
# List current assemblies which have never been imported
python bin/registry_tracker.py \
  --host mysql.registry.example.org \
  --port 4527 \
  --user registryuser \
  --password registry-secret \
  --database refseq_assembly \
  pending

# Reserve an assembly before starting the pipeline; capture the printed id
status_id=$(python bin/registry_tracker.py \
  --host mysql.registry.example.org \
  --port 4527 \
  --user registryuser \
  --password registry-secret \
  --database refseq_assembly \
  start \
  --accession GCF_000001405.40 \
  --genebuilder automated)

# Call this only after the import has completed successfully
python bin/registry_tracker.py \
  --host mysql.registry.example.org \
  --port 4527 \
  --user registryuser \
  --password registry-secret \
  --database refseq_assembly \
  complete \
  --import-status-id "$status_id" \
  --genebuild-version EXT01 \
  --last-genebuild-update 2026-01
```

`start` refuses to create a row when any import history already exists for the
current assembly. The Nextflow pipeline uses `automated` for `genebuilder` by
default. Direct calls to `registry_tracker.py` must provide
`--genebuilder NAME`.

By default, an existing import-status record stops the pipeline. Set
`--skip_start_check true` (or the bare `--skip_start_check` flag) to reuse that
record, but only when its
`genebuilder` matches the requested value. Reuse resets the record to
`in_progress` and updates its start and status-update dates.

Static core metakeys are configured in
[`config/static_metakeys.json`](config/static_metakeys.json). Entries under
the `refseq` key are inserted into the core `meta` table by
`ADD_STATIC_METAKEYS`. Existing metakeys are preserved because the insert uses
`INSERT IGNORE`.

## Next steps
Run the statistics pipeline on the output CSV to generate statistics for handover.

Create a params file with all the analyses you wish to run. 
For handover, you want to run all analyses.

```bash
{
  "csvFile": "/path/to/statistics_input.csv",
  "outdir": "stats_out",
  "run_busco_core": true,
  "busco_mode": "both",
  "apply_busco_metakeys": true,
  "run_ensembl_stats": true,
  "apply_ensembl_stats": true,
  "run_ensembl_beta_metakeys": true,
  "apply_ensembl_beta_metakeys": true,
  "run_pepstats": true,
  "team": "genebuild",
  "host": "",
  "port": 4532,
  "user": "",
  "user_r": "",
  "password": "",
  "enscode": "/path/to/enscode"
}
```

Once you have a params file, you can run the statistics pipeline.

```bash
nextflow run pipelines/statistics/main.nf \
  -params-file /path/to/stats_params.json
```

## Run production sync as normal
