# Public channel contracts

All public tuples retain the first `meta` map. Transformed maps use
`meta + [new_key: value]`; shared metadata is never mutated in place.

## Named subworkflow boundaries

| Subworkflow | Public inputs | Public outputs |
|---|---|---|
| `PREPARE_LONG_READS` | manifest, input-run requests, preparation policy | `reads`, normalized manifest/checksums/status, versions |
| `ALIGN_LONG_READS` | prepared reads, reference, alignment policy | `alignment`, QC, features, alignment manifests, reference dictionary/identity, versions |
| `PREPARE_CONTIG_SHARDS` | indexed alignments, shard policy | workload, plans, shard tuples, summaries/status, versions |
| `AUDIT_NATIVE_SHARDS` | native model tuples, empty-plan records, shard plans, backend list | audited native tuples, backend-plan audit versions |
| `MERGE_LONG_READ_MODELS` | audited native tuples, scope policy, merge arguments | normalized-input GTF tuples with file lists/checksums, versions |
| `NORMALIZE_VALIDATE_MODELS` | merged GTF tuples, reference dictionary | validated normalized-model tuples, reports, versions |
| `COLLECT_LONG_READ_EVIDENCE_PRODUCTS` | validated models, normalized manifest, alignment manifests, reference identity | verified evidence package, versions |
| `RUN_DIAMOND_QC` | validated normalized models and optional protein database | backend/scope-separated ORF, DIAMOND, and comparison products |
| `FINALIZE_LONG_READ_PIPELINE` | normalized inputs, shard plans, model artifacts, published-model mapping, merge provenance, reference identity, stage/evidence status | run manifest, final checksums, checksum verification, software versions |

The deterministic fixture exercises these boundaries through the all-backend
stub matrix and the real host-tool tmerge fixture; Python tests cover the
serialization and validation helpers independently.

* Input run: `tuple(meta, fastq)` or `tuple(meta, bam, bai)`; `meta` contains `id`, `cohort_id`, `run_id`,
  `technology`, `assembly_id`, `read_set_id`, `sample_id`, `library_id`,
  `read_representation`, and `source_uri`. `raw_subreads` is accepted only
  when `ccs_mode=pbccs`; it must carry a PBI sidecar and is replaced by the
  generated CCS BAM before alignment.
  When `preparation_stages` includes `acquire`, the preparation request is
  `tuple(meta, source_uri)` and acquisition emits the canonical artifact tuple
  after a non-empty download/copy plus a checksum-bearing manifest row.
  `bam_to_fastq` emits canonical FASTQ before FASTQ validation and staging.
* Alignment: `tuple(meta, bam, bai)`. QC, features, versions, and an alignment
  manifest are separate named outputs. The manifest records reference and BAM
  checksums plus minimap2 preset and secondary-alignment policy.
* Workload: `tuple(meta, bam, bai, workload_tsv)`.
* Shard: `tuple(meta, shard_meta, shard_bam, shard_bai)`; `shard_meta` contains
  `shard_id`, `contig`, `reference_length`, `mapped_reads`, `resource_class`,
  `source_run_id`, and `shard_status`.
  `VALIDATE_SHARD_COVERAGE` checks that emitted shard mapped-read totals equal
  the source BAM; zero-read plans are recorded without launching tool tasks.
* Native model: `tuple(meta, backend_meta, native_product, status_report)`.
* Normalized model: `tuple(meta, backend, scope, bed12, model_manifest,
  validation_report)`.

The shard plan is authoritative for planned, emitted, empty, skipped, and
failed shards. Reference dictionary compatibility is checked before the plan
is created. Native formats are backend-specific and cannot enter comparison or
evidence collection without conversion to validated BED12. The backend audit
compares eligible plan IDs with native audit records and records skipped,
successful, empty, and failed/missing states.

Final model checksum rows use `path` for the eventual published product path and
`source_path` for the staged file verified during finalization, because Nextflow's
`publishDir` copy can complete asynchronously after the finalizer task.
