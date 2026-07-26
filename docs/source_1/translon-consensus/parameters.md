# Translon-Consensus Parameters

Automatically generated from the pipeline `nextflow_schema.json`.

## Parameters

Parameter schema for the translon consensus pipeline

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `bed_results_dir` | string |  | yes | Path to directory containing BED12 results organized by tool (e.g., results/*/sample.bed12). BED files must be in proper BED12 format with one entry per ORF. Use conversion scripts in scripts/ directory to convert tool-specific formats. |
| `gencode_gtf` | string |  | no | Path to gencode GTF annotation file for adding RNA biotype and CDS context. Can be gzipped (.gtf.gz). If not provided, features will not have biotype/CDS annotations. |
| `gencode_fasta` | string |  | no | Path to gencode FASTA. |
| `gencode_fasta_fai` | string |  | no | Path to gencode FASTA FAI. |
| `ucsc_session_url` | string | https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38 | no | UCSC Genome Browser session URL (without position) for generating genome browser links in reports |
| `outdir` | string | ./results | no | Output directory for pipeline results |
| `tracedir` | string |  | no | Directory for pipeline execution reports |
| `ftp_dir` | string |  | no | Optional FTP directory path for moving the HTML report. If not provided, report will only be saved locally. |
