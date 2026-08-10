# ensembl_import Handover

## Current State

The pipeline lives at `pipelines/ensembl_import`. It was previously under the misspelled directory `pipelines/ensebl_import`; use the corrected path from now on.

The RefSeq import path is wired and stub-tested with Nextflow 25.04.6:

```bash
env NXF_TEMP=/tmp \
  /hps/software/spack/opt/spack/linux-cascadelake/nextflow-25.04.6-og6z2234zw6ecuuv44fyaaynvqwkxlcu/bin/nextflow \
  run main.nf \
  -stub \
  --input_csv /tmp/ensembl_import_stub/input.csv \
  --server_settings /tmp/ensembl_import_stub/server_settings.json \
  --ensembl_genes_repo /hps/software/users/ensembl/genebuild/lazar/modenv/stats_pipe/ensembl-genes-nf \
  --outdir /tmp/ensembl_import_stub/results
```

Slurm is the default executor for this pipeline. Use `-profile local` only for local development or stub checks.

The registry update path is not implemented yet. The registry modules are placeholders and should not be included in the runnable path until they have valid process names, valid input declarations, outputs, versions, and stubs.

## Required Nextflow Version

Use Nextflow 25.x. The available module on this system is `nextflow/25.04.6`.

The root `nextflow.config` pins `manifest.nextflowVersion = '25.04.6'`, and `pipelines/ensembl_import/main.nf` has a runtime guard that rejects non-25.x versions.

The default `nextflow` currently resolves to 26.04.0 on this host, so do not rely on PATH unless the 25.x module is loaded first.

## Execution

`pipelines/ensembl_import/nextflow.config` sets `process.executor = 'slurm'`, so normal runs submit process jobs to Slurm.

Use the default Slurm execution for real imports:

```bash
env NXF_TEMP=/tmp \
  /hps/software/spack/opt/spack/linux-cascadelake/nextflow-25.04.6-og6z2234zw6ecuuv44fyaaynvqwkxlcu/bin/nextflow \
  run /hps/software/users/ensembl/genebuild/lazar/modenv/stats_pipe/ensembl-genes-nf/pipelines/ensembl_import/main.nf \
  --input_csv input.csv \
  --ensembl_genes_repo /hps/software/users/ensembl/genebuild/lazar/modenv/stats_pipe/ensembl-genes \
  --server_settings server_settings.json \
  --outdir /hps/nobackup/flicek/ensembl/genebuild/lazar/ensembl_load_test/pipeline
```

Use local execution only when requested explicitly:

```bash
env NXF_TEMP=/tmp \
  /hps/software/spack/opt/spack/linux-cascadelake/nextflow-25.04.6-og6z2234zw6ecuuv44fyaaynvqwkxlcu/bin/nextflow \
  run main.nf \
  -profile local \
  -stub \
  --input_csv /tmp/ensembl_import_stub/input.csv \
  --server_settings /tmp/ensembl_import_stub/server_settings.json \
  --ensembl_genes_repo /hps/software/users/ensembl/genebuild/lazar/modenv/stats_pipe/ensembl-genes-nf \
  --outdir /tmp/ensembl_import_stub/results
```

## RefSeq Import Flow

`main.nf` should stay thin:

1. Validate parameters with `nf-schema` plus local file existence checks.
2. Read `--input_csv` with columns `gcf` and `species`.
3. Build one metadata map per assembly: `[id: row.gcf, species: row.species]`.
4. Read `--server_settings` JSON once and pass database settings as a value channel.
5. Call `IMPORT_REFSEQ`.

`IMPORT_REFSEQ` coordinates these modules:

1. `FETCH_REFSEQ`: runs the Ensembl `gff_cli.py refseq run` command and emits GFF3, FASTA, and assembly report files.
2. `LOAD_REFSEQ`: runs `gff_cli.py create-core` and emits a done marker so downstream steps depend on successful loading.
3. `GET_SAMPLE_GENE`: runs `bin/get_sample_gene.py` against the loaded core database and inserts the six `sample.*` keys into the core `meta` table.

`get_sample_gene.py` is Python 3.10-compatible and follows the Ensembl `HiveAddPlaceholderLocation` selection rules: it searches the ten longest sequence regions for a supported protein-coding transcript (`hcoverage >= 99`, `percent_id >= 75`), or otherwise chooses the protein-coding gene with the most exons.

The script can also be run directly after a core database has been loaded:

```bash
python pipelines/ensembl_import/bin/get_sample_gene.py \
  --db-name species_core_114_1 \
  --host mysql.example.org \
  --port 3306 \
  --user ensrw \
  --password 'PASSWORD'
```

The database credentials are passed from `--server_settings` when the Nextflow process runs.

## Repository Pipeline Rules

Follow the rules from `docs/template`:

- Use DSL2 in every runnable workflow/subworkflow file.
- Keep `main.nf` as orchestration only; implementation logic belongs in modules.
- Structure pipelines as `main.nf -> subworkflows -> modules`.
- Name subworkflows and processes in uppercase; name files in lowercase with underscores.
- Use `../modules/...` includes from pipeline subworkflows.
- Use named `emit` outputs for every module and subworkflow output.
- Emit `versions.yml` from every module.
- Add `stub` blocks for every module so `nextflow run main.nf -stub` validates channel wiring quickly.
- Use resource labels such as `process_light`, `process_low`, `process_medium`, and `process_high_memory`.
- Keep tool arguments configurable with `task.ext.args` and optional execution with `task.ext.when`.
- Preserve metadata maps through sample-like data channels.
- Do not use single-element `tuple val(meta)` inputs; use `val meta` for metadata-only channels.
- Use tuples when passing metadata plus files or markers, for example `tuple val(meta), path(file)`.
- Validate required CSV columns before the pipeline launches work; this pipeline expects non-empty `gcf` and `species`.
- Do not hardcode local paths in modules; use params such as `params.ensembl_genes_repo` and config overrides.
- Publish outputs under pipeline-specific subdirectories inside `params.outdir`.

## Inputs

`--input_csv` must contain:

```csv
gcf,species
GCF_000001405.40,Homo sapiens
```

`--server_settings` must be JSON with:

```json
{
  "db_host": "localhost",
  "db_port": 3306,
  "db_user": "ensro",
  "db_password": "",
  "db_read_user": "ensro"
}
```

`--ensembl_genes_repo` must point to a local checkout containing:

- `src/python/ensembl/genes/ensembl_loading/gff_cli.py`
- `src/python/ensembl/genes/metadata/core_metadata.py`

## Follow-Up Work

- Implement the registry update subworkflow after the registry API and required database actions are confirmed.
- Replace placeholder version values for Ensembl Python scripts with real version extraction if those scripts expose a version command.
- Run a non-stub RefSeq import against a test database once `ensembl_genes_repo`, RefSeq network access, and MySQL credentials are available.
