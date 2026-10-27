# QC Pipeline

This pipeline runs annotation QC from a genome FASTA and a GFF3 annotation. It
can run AGAT metrics, InterProScan protein coverage, Diamond protein
validation, or all three.

## Requirements

- Nextflow 26.04.6 or later (strict syntax)
- Singularity or Apptainer
- The InterProScan 5.78-109.0 data directory for `interproscan` and `combined`
- A QC parser container exposing `annotation-qc parse-agat` and
  `annotation-qc parse-interpro` when using AGAT or InterProScan
- The pipeline supplies its Ensembl AGAT feature definitions internally
- Registry credentials configured outside the pipeline if the GitLab registry
  requires authentication

The parser image is assigned centrally to the `qc_parser` process label and
currently uses the EBI Docker registry `:latest` image. Once the pipeline is
versioned, change that single config reference to an immutable container tag
based on the corresponding commit SHA.

To select the Nextflow version explicitly, set `NXF_VER` before the command:

```bash
NXF_VER=26.04.6 nextflow run pipelines/qc/main.nf --help
```

## Input

The manifest columns depend on the selected mode. The pipeline validates the
manifest against a mode-specific schema before starting any processes.

For AGAT, only the GFF3 is required:

```csv
sample,gff3
sample_1,/path/to/sample_1.gff3
```

For InterProScan, either a precomputed protein FASTA or both the genome FASTA
and GFF3 may be supplied. The latter allows proteins to be derived with
`gffread`:

```csv
sample,protein_fasta
sample_1,/path/to/sample_1.faa
```

For Diamond and combined mode, the genome FASTA and GFF3 are required because
Diamond translation validation uses both. A precomputed protein FASTA is
optional in those modes:

```csv
sample,genome_fasta,gff3,protein_fasta
sample_1,/path/to/sample_1.fa,/path/to/sample_1.gff3,/path/to/sample_1.faa
```

InterProScan and Diamond use the supplied protein FASTA when present and
otherwise derive one with `gffread` where the required inputs are available.

Diamond requires either a prebuilt database via `--diamond_reference_db` or a
reference protein FASTA via `--diamond_reference_proteins`. The prebuilt
database takes precedence and avoids rebuilding it. The FASTA is used only as
the fallback source for `diamond makedb`.

## Parameters

| Parameter | Required | Description |
| --- | --- | --- |
| `--input_csv` | yes | Mode-specific CSV manifest; see the Input section |
| `--mode` | yes | `agat`, `interproscan`, `diamond`, or `combined` |
| `--data_file_path` | InterProScan modes | InterProScan data directory; supply this through a site/user config |
| `--diamond_reference_db` | Diamond modes | Existing prebuilt Diamond database; takes precedence over the FASTA |
| `--diamond_reference_proteins` | Diamond modes | Reference protein FASTA used to build a database when no prebuilt database is supplied |
| `--database` | no | InterProScan database, default `Pfam` |
| `--outdir` | no | Output directory, default `./results` |

## Running

```bash
nextflow run pipelines/qc/main.nf \
  -profile slurm \
  -c ~/.config/ensembl-genes-nf/ebi.config \
  --input_csv /path/to/annotation_samples.csv \
  --mode combined \
  --diamond_reference_db /path/to/reference.dmnd \
  --outdir results
```

The InterProScan data directory is site-specific and is intentionally not
hard-coded in the repository. Create a user or site configuration containing,
for example:

```groovy
params.data_file_path = '/nfs/production/flicek/ensembl/shared_data/interproscan-5.78-109.0/data'
```

Pass that file with `-c`. The parameter schema requires the value for
`interproscan` and `combined` modes and validates that the directory exists.

Do not put registry usernames, passwords, or tokens in the sample sheet or
committed Nextflow configuration. Configure Apptainer or Singularity
authentication on the execution host before launching the run.

## Registry Authentication

The QC parser image is private in the EBI Docker registry. Authenticate
Apptainer once on the execution host using an EBI/GitLab account or token with
permission to read the registry:

Enter the token interactively and pipe it to the login command. The token is
never placed on the command line or stored in the repository:

```bash
export EBI_REGISTRY_USER="${EBI_REGISTRY_USER:-$USER}"
read -r -s EBI_REGISTRY_TOKEN
printf '\n'
printf '%s\n' "$EBI_REGISTRY_TOKEN" | singularity registry login \
  --username "$EBI_REGISTRY_USER" \
  --password-stdin \
  docker://dockerhub.ebi.ac.uk
unset EBI_REGISTRY_TOKEN
```

Verify the credentials and image access before starting Nextflow:

```bash
singularity registry list
singularity pull docker://dockerhub.ebi.ac.uk/ensembl_genebuild/ensembl-genes-containers/ensembl-genes-qc-dev:latest
```

For CI or other non-interactive runs, provide masked secrets as
`APPTAINER_DOCKER_USERNAME` and `APPTAINER_DOCKER_PASSWORD`. Do not pass
registry credentials as Nextflow parameters or commit them to configuration.

## Outputs

Published outputs are written to:

```text
<outdir>/qc/agat/
<outdir>/qc/interpro/
<outdir>/qc/diamond/
```

When Diamond builds a database from FASTA, the reusable database is also
published as `<outdir>/qc/diamond/reference.dmnd` and can be supplied to a
later run with `--diamond_reference_db`.

## Separation Of Concerns

- `nextflow_schema.json` validates parameter types, enums, and required
  mode-specific parameters.
- `assets/*_samplesheet.json` validates the mode-specific annotation manifest.
- `main.nf` loads validated inputs and starts the workflow.
- `workflows/qc.nf` decides which analyses run and shares derived proteins.
- Subworkflows connect reusable analysis chains.
- Modules run one tool or parser and do not decide which QC mode the user chose.

In short: **the schema validates configuration; workflow code resolves that
configuration into behaviour.** This keeps reusable subworkflows independent
of pipeline-specific parameter names and prevents one universal sample sheet
from accumulating unrelated future inputs.
