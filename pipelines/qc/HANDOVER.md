# QC Pipeline Handover

## Current State

The QC pipeline runs AGAT metrics when `--run_agat_metrics` is enabled and can optionally run InterProScan with `--run_interproscan`.

The InterPro branch was repaired to match the repository module pattern:

- input rows are split into separate GFF3 and protein channels
- the InterPro process now has valid syntax, output naming, `versions.yml`, and a `stub` block
- the parser process now accepts the parser path correctly and also emits `versions.yml`
- the subworkflow now wires the process outputs together correctly

The AGAT side was brought up to the same module standard so the pipeline is consistent end to end.

## What Was Fixed

- `run_interpro.nf`
  - fixed the broken `containerOptions` block
  - pinned the module to Docker image `interpro/interproscan:5.78-109.0`, executed through Singularity
  - uses the full `/opt/interproscan/interproscan.sh` path required by Singularity
  - fixed the command line continuation bugs in the `script` block
  - fixed output naming so it uses the sample/stem instead of stringifying the metadata map
  - added `versions.yml`, `when`, `task.ext.args`, and `stub`

- `interpro_stats.nf`
  - fixed the missing `interpro_parser` input
  - fixed the undefined parser variable
  - passes the protein FASTA from the input CSV to `interpro.py --query_protein`
  - publishes the parsed result as `<sample>_hits.tsv` directly under `qc/interpro/`
  - added `versions.yml`, `when`, `task.ext.args`, and `stub`

- `subworkflows/interpro/interproscan.nf`
  - fixed the channel wiring between the runner and parser
  - added combined version emission

- `main.nf`
  - split the sample sheet into `gff_ch` and `protein_ch`
  - added validation that rows contain the file type required by the enabled branch
  - passed `interpro_parser` through to the InterPro subworkflow

- `nextflow.config`
  - restored `run_agat_metrics = false` so InterPro-only runs do not require AGAT inputs and feature-level config

- `README.md` and `nextflow_schema.json`
  - updated the documented input CSV columns and InterPro parameters

## Execution Notes

- Use Nextflow 25.x for this pipeline.
- Use InterProScan 5.78-109.0 for the QC InterPro branch. InterProScan 6 is a separate Nextflow workflow and is not compatible with this module's `interproscan.sh` interface.
- Do not rely on the host `nextflow` PATH default if it resolves to 26.x.
- Input CSV format is now:

```csv
sample,gff3,protein
sample_1,/path/to/sample_1.gff3,/path/to/sample_1.faa
```

- `--input_csv` is required.
- `--ensembl_genes_repo` is required.
- `--run_agat_metrics` defaults to false.
- `--run_interproscan` defaults to false.
- `--database` defaults to `Pfam`.
- `--data_file_path` should point at the InterProScan data directory when the InterPro branch is enabled.
- Run with `-profile slurm` to submit QC tasks to Slurm and execute the Docker image through Singularity.
- InterProScan temporary files are written under `<outdir>/qc/interpro/tmp/<sample>/`, rather than the Nextflow work directory.

## Follow-Up Work

- Run a stubbed Nextflow check against the QC pipeline once the desired Nextflow binary is available on the host.
- If the upstream `ensembl-genes` parser path changes, update the fallback path in `interpro_stats.nf` and the README.
- Consider adding a small pipeline-level integration test fixture with one GFF3 and one protein FASTA so both branches can be exercised together.
