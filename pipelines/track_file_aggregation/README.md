# Track-file aggregation

This DSL2 pipeline is an optional downstream stage for
`pipelines/track_file_preparation`. It consumes the preparation pipeline's
`track_manifest.tsv` plus an explicit aggregation manifest. It never discovers
entities or infers sample/tissue membership from filenames or directories.

The pipeline supports two explicit aggregation modes:

- `raw`: sum source coverage values and splice-junction read counts;
- `equal_sample`: sum after assigning each contributing sample a weight of
  `1/N`. Coverage values are scaled before BigWig conversion and junction
  counts are rounded to integer support for the BigBed schema.

Gene-model BigBeds are intentionally not aggregated. Combining annotations
requires a separate merge and deduplication policy.

## Inputs

The aggregation manifest has one row per source entity and these required
columns:

```text
aggregate_id  aggregate_type  source_entity_id  source_entity_type  aggregation_method
```

Optional columns are `gca_accession`, `aggregate_label`, and `metadata_json`.
`aggregation_method` is `raw` or `equal_sample`; the latter requires
`source_entity_type=sample`. Aggregate types are `sample`, `tissue`, and
`global`. Tissue membership is explicit: callers must assign missing tissue
metadata to an aggregate such as `uncharacterized`; the workflow does not
infer tissue from `organism_part` or other source attributes. One assembly is
required per invocation. Use two invocations for runs-to-samples followed by
samples-to-tissues.

Coverage aggregation operates on existing normalized BigWigs; it does not pool
reads or recalculate coverage from BAMs. The preparation manifest must contain
machine-readable JSON in `normalization_parameters` so incompatible bin sizes
or normalization methods are rejected. Aggregate manifests preserve this
metadata so sample-level outputs can be consumed by a later tissue-level run.

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
