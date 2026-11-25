# Supporting evidence aggregation

This DSL2 pipeline creates a per-assembly evidence bundle for a custom QC
reader. It consumes the run-level `track_manifest.tsv` produced by
`track_file_preparation` and keeps one source entity as one set of columns.
It does not infer sample, tissue, or cohort membership and it does not sum
support across runs.

The invocation is assembly-scoped because BigWig and BigBed coordinates must
share one reference sequence dictionary. Run IDs are retained as the column
labels. Pass `--species_name` when the assembly accession is not the preferred
QC label. The pipeline validates source checksums, source track schemas, and
coverage normalization metadata before any matrix task starts.

Coverage produces one deepTools matrix and raw-count table per strand. The
matrix has one column per source BigWig. Splice junctions are converted from
the preparation BigBeds and written as a sparse wide TSV with two columns per
source entity, `<run>.unique_reads` and `<run>.multi_reads`. Missing junctions
are represented by zeroes while a missing source track is a validation error.

Outputs are published under:

```text
<outdir>/plan/
<outdir>/coverage/
<outdir>/splice_junctions/
<outdir>/reports/
```

The report directory contains `supporting_evidence_manifest.tsv`,
`supporting_evidence_status.tsv`, and
`supporting_evidence_metadata.json`. These files provide the QC reader with
the assembly, source entities, output paths, and checksums.

Example:

```bash
nextflow run pipelines/supporting_evidence_aggregation/main.nf \
  -profile test \
  --track_manifest pipelines/track_file_aggregation/test/track_manifest.tsv \
  --manifest_base_dir pipelines/track_file_aggregation/test \
  --outdir /tmp/supporting-evidence-aggregation-stub
```

For a production run, use a matching assembly-specific track manifest and run
with the repository's HPC Singularity configuration. The source manifest is
the contract between track preparation and this pipeline.
