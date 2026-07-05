# Repeat Annotation Pipeline

The repeat annotation pipeline fetches genome assemblies, prepares RepeatModeler libraries, and runs selected repeat annotation tools. It can reuse a pre-computed RepeatModeler library from Ensembl FTP or generate a new library when one is missing.

## Entry Point

```bash
nextflow run pipelines/repeat/main.nf \
  -c pipelines/repeat/nextflow.config \
  --csvFile repeat_samples.csv \
  --outdir results/repeat \
  --generate_lib \
  --run_repeatmasker
```

Enable additional tools with `--run_red`, `--run_dust`, and `--run_trf`.

## Input CSV

The CSV is validated by `pipelines/repeat/assets/schema_input.json`.

| Column | Required | Description |
| --- | --- | --- |
| `species_name` | yes | Species name without spaces. Used in Ensembl repeat FTP paths. |
| `gca` | yes | Genome assembly accession such as `GCA_000001405.29`. |
| `genome_file` | no | Optional local genome FASTA path. Leave blank to fetch by accession. |

Example:

```text
species_name,gca,genome_file
homo_sapiens,GCA_000001405.29,
```

## Workflow Stages

| Stage | Module | Purpose |
| --- | --- | --- |
| Fetch genome | `FETCH_GENOME` | Download the assembly or use a supplied genome file. |
| Check repeat library | `FETCH_REPEAT_MODEL` | Look for an existing RepeatModeler library on FTP. |
| Generate library | `GENERATE_REPEATMODELER_LIBRARY` | Build a de novo RepeatModeler library when needed. |
| Download library | `CHECK_AND_DOWNLOAD_RMLIBRARY` | Download an existing library when available. |
| Upload library | `UPLOAD_INTO_FTP` | Publish newly generated libraries to the configured FTP target. |
| Annotate repeats | `RUN_REPEATMASKER`, `RUN_RED`, `RUN_DUST`, `RUN_TRF` | Run selected repeat annotation tools. |

## Output Layout

Outputs are published per assembly accession:

```text
results/repeat/
└── <GCA>/
    ├── library/
    ├── rm_library/
    ├── repeatmasker/
    ├── red/
    ├── dust/
    └── trf/
```

Tool outputs currently focus on GTF files and `versions.yml` provenance files. The shared execution reports are written to `results/pipeline_info` unless `--tracedir` is overridden.

## Parameters

```{toctree}
:maxdepth: 1

../generated/repeat-parameters
```
