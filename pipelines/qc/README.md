# QC Pipeline

This pipeline runs annotation quality-control metrics for one or more GFF3 files.

The implemented QC steps are:

- AGAT statistics generation followed by conversion of the AGAT text report into a `genebuild` CSV using the `parse_agat.py` parser from the `ensembl-genes` repository.
- Pairwise annotation comparison for two GFF3 annotation sources on the same assembly.

The pairwise branch is an initial port of the HPRC Ensembl/CAT comparison scripts. Internally, some output columns still use the legacy `ensembl`/`cat` names: source A maps to the legacy Ensembl side, and source B maps to the legacy CAT side.

## What It Does

For each row in an AGAT input CSV, the pipeline:

1. runs `agat_sp_statistics.pl` on the input GFF3 file
2. writes the AGAT text report
3. parses that report into a `*_genebuild.csv` metrics file

For each row in a pairwise input CSV, the pipeline:

1. stages the two GFF3 inputs from local paths or URLs
2. computes gene-level overlaps and reciprocal best-hit pairs
3. optionally runs transcript concordance, CDS integrity, gene presence, multi-mapping, feature count, transcript count, and bidirectional `gffcompare` metrics

## Requirements

- Nextflow with DSL2 enabled
- Singularity or Apptainer available on the execution host
- For AGAT metrics, a local checkout of `ensembl-genes` containing:
  - `src/python/ensembl/genes/annotation-qc/parsers/parse_agat.py`
  - `src/python/ensembl/genes/annotation-qc/metrics/config/feature_levels.yaml`

## Inputs

### AGAT Sample Sheet

The AGAT input is `--gff_csv`, a CSV file with a header and at least these columns:

```csv
sample,gff3
sample_1,/path/to/sample_1.gff3
sample_2,/path/to/sample_2.gff3
```

- `sample`: sample identifier used for logging and task tags
- `gff3`: absolute or relative path to the GFF3 file

### Pairwise Sample Sheet

Enable the pairwise branch with `--run_pairwise_annotation_comparison true` and pass `--pairwise_csv`.

The pairwise CSV supports local paths and `http`, `https`, or `ftp` URLs interchangeably:

```csv
sample,source_a,source_b,gff_a,gff_b,assembly_report
sample_1,ensembl,refseq,/path/to/ensembl.gff3,https://example.org/refseq.gff3.gz,/path/to/assembly_report.txt
```

- `sample`: sample identifier used for logging and output paths
- `source_a`, `source_b`: source labels; optional, default to `source_a` and `source_b`
- `gff_a`, `gff_b`: local path or URL for each annotation GFF3
- `assembly_report`: optional local path or URL reserved for contig-normalization workflows

URL inputs are checked against `--pairwise_cache_dir` by basename before downloading. Once staged, downstream modules consume fixed local filenames, so local paths and URLs follow the same workflow.

## Parameters

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--gff_csv` | when AGAT enabled | none | CSV describing the samples and GFF3 paths |
| `--ensembl_genes_repo` | when AGAT enabled | none | Path to a local `ensembl-genes` checkout |
| `--feature_levels` | no | derived from `ensembl_genes_repo` | Path to `feature_levels.yaml` |
| `--agat_parser` | no | derived from `ensembl_genes_repo` | Path to `parse_agat.py` |
| `--run_agat_metrics` | no | `true` | Enable the AGAT metrics branch |
| `--run_pairwise_annotation_comparison` | no | `false` | Enable the pairwise annotation comparison branch |
| `--pairwise_csv` | when pairwise enabled | none | CSV describing source A/source B annotation pairs |
| `--pairwise_cache_dir` | no | `<outdir>/cache/pairwise_annotation` | Cache directory for URL inputs |
| `--pairwise_contig_normalization` | no | `none` | Legacy gene-pairing contig normalization mode: `none`, `basic`, or `cat_hash` |
| `--pairwise_ensg_lookup` | no | none | Optional lookup TSV for the legacy gene-presence parser |
| `--run_pairwise_gffcompare` | no | `false` | Enable bidirectional `gffcompare` metrics |
| `--outdir` | no | `./results` | Output directory |

## Running

Run AGAT metrics:

```bash
nextflow run pipelines/qc/main.nf \
  --gff_csv /path/to/gff_samples.csv \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --outdir results \
  -process.executor slurm
```

Run only pairwise annotation comparison:

```bash
nextflow run pipelines/qc/main.nf \
  --run_agat_metrics false \
  --run_pairwise_annotation_comparison true \
  --pairwise_csv /path/to/pairwise_samples.csv \
  --outdir results
```

## Outputs

AGAT outputs are written to:

```text
<outdir>/qc/agat/
```

Per sample, AGAT produces:

- `<sample>_agat_stats.txt`: raw AGAT statistics output
- `<sample>_agat_stats_genebuild.csv`: parsed metrics in CSV format

Pairwise outputs are written under:

```text
<outdir>/qc/pairwise_annotation/<sample>/
```

The first-pass port can produce:

- `gene_pairs/*.gene_pairs_all.tsv`
- `gene_pairs/*.gene_pairs_rbh.tsv`
- `gene_pairs/*.assembly_summary.tsv`
- `transcript_concordance/*_transcript_concordance.tsv`
- `coding_integrity/*_coding_integrity.tsv`
- `gene_presence/*_gene_presence.tsv`
- `multi_mapping/*_multi_mapping.tsv`
- `feature_metrics/*.tsv`
- `transcript_counts/*_gene_transcript_counts.tsv`
- `gffcompare/*` when `--run_pairwise_gffcompare true`

## Validation

The workflow stops early if required inputs for the enabled branches are missing or invalid.

## Implementation Notes

- The workflow entrypoint is `pipelines/qc/main.nf`.
- AGAT statistics are implemented in `modules/agat/run_agat_stats.nf`.
- AGAT report parsing is implemented in `modules/agat/parse_agat.nf`.
- The AGAT subworkflow is defined in `subworkflows/agat/agat_stats.nf`.
- Pairwise annotation comparison is defined in `subworkflows/pairwise_annotation/pairwise_annotation_comparison.nf`.
- Pairwise modules are under `modules/pairwise_annotation/`.

## Current Scope

This README reflects the current implementation in this branch. The pairwise comparison path is a mechanical port and still needs a schema cleanup pass to replace legacy Ensembl/CAT column names with generic source A/source B naming.
