# Qc Parameters

Automatically generated from the pipeline `nextflow_schema.json`.

## Parameters

Parameter schema for the QC pipeline

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `gff_csv` | string |  | yes | Path to a CSV file with at least sample and gff3 columns |
| `ensembl_genes_repo` | string |  | yes | Path to a local ensembl-genes checkout |
| `feature_levels` | string |  | no | Optional explicit path to feature_levels.yaml |
| `agat_parser` | string |  | no | Optional explicit path to parse_agat.py |
| `run_agat_metrics` | boolean | True | no | Enable the AGAT metrics branch |
| `outdir` | string | ./results | no | Output directory |
| `tracedir` | string |  | no | Directory for Nextflow execution reports and trace files |
