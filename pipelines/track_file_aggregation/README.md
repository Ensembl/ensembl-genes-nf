# Track-file aggregation

This DSL2 pipeline is an optional downstream stage for
`pipelines/track_file_preparation`. It consumes the preparation pipeline's
`track_manifest.tsv` plus an explicit aggregation manifest. It never discovers
entities or infers sample/tissue membership from filenames or directories.

The first implementation aggregates:

- forward and reverse coverage BigWigs by equal-weight mean;
- splice-junction BigBeds by summing unique and multimapping read counts.

Gene-model BigBeds are intentionally not aggregated. Combining annotations
requires a separate merge and deduplication policy.

## Inputs

The aggregation manifest has one row per source entity and these required
columns:

```text
aggregate_id  aggregate_type  source_entity_id  source_entity_type
```

Optional columns are `gca_accession`, `aggregate_label`, `coverage_method`,
`splice_method`, and `metadata_json`. `coverage_method=mean` and
`splice_method=sum_counts` are currently supported. One assembly is required
per invocation. Use two invocations for runs-to-samples followed by
samples-to-tissues.

Coverage aggregation averages existing normalized BigWigs; it does not pool
reads or recalculate coverage from BAMs. The preparation manifest must contain
machine-readable JSON in `normalization_parameters` so incompatible bin sizes
or normalization methods are rejected.

## Example

```bash
nextflow run pipelines/track_file_aggregation/main.nf \
  -profile test \
  --track_manifest pipelines/track_file_aggregation/test/track_manifest.tsv \
  --aggregation_manifest pipelines/track_file_aggregation/test/aggregation_manifest.tsv \
  --chrom_sizes pipelines/track_file_aggregation/test/chrom.sizes \
  --outdir /tmp/track-file-aggregation-stub
```

Outputs include `track_manifest.tsv`, `aggregate_status.tsv`,
`exceptions.tsv`, checksums, provenance JSON, and the aggregate track files
under the configured output directory.
