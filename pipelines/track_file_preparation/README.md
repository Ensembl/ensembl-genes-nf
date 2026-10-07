# Track-file preparation

This DSL2 pipeline treats the supplied TSV as the authoritative entity
inventory. It does not scan production directories or infer run/sample/merged
relationships. The per-row `track_types` value is used unless the pipeline
parameter supplies an explicit override; `all` means coverage, gene-model,
and splice-junction tracks. Junction derivation from BAM is opt-in with
`--derive_junctions true`.

Run the wiring test with:

```bash
nextflow run pipelines/track_file_preparation/main.nf -stub-run -profile test \
  --input_manifest pipelines/track_file_preparation/test/entity_manifest.tsv \
  --chrom_sizes pipelines/track_file_preparation/test/test.chrom.sizes \
  --assembly_release test-release --outdir /tmp/track-file-preparation-stub
```

The manifest must contain the required columns from
`TRACK_FILE_PREPARATION_SPEC.md`. Relative paths are resolved against the
manifest directory or `--manifest_base_dir`. Output records retain the
original entity ID while filenames use a collision-checked sanitized ID.

For BAM-dependent tracks, the pipeline uses a manifest-supplied BAM index or
an adjacent `.bai`/`.csi` when available. If neither exists, it creates a
coordinate-sorted BAM `.bai` with `samtools index` before running coverage or
BAM-derived junction extraction. Index creation is a declared workflow stage,
not a prerequisite that the manifest producer must perform.
