# Long-read reconstruction pipeline

This is the strict DSL2 replacement for `pipelines/long_read_tama/`. It accepts
explicit FASTQ, indexed BAM, approved-manifest, and retrospective inputs and
publishes normalized BED12 products under run- and cohort-namespaced paths.

The public contracts are documented in `docs/CHANNEL_CONTRACTS.md`. Required
stages fail the invocation; a selected backend is never converted into an
optional best-effort task. Native backend products are retained only when
`--publish_native_intermediates true` is selected. The legacy pipeline is not a
runtime dependency.

Preparation extensions are explicit through `preparation_stages`: `discover`,
`resolve_metadata`, `acquire`, and `bam_to_fastq`. Acquisition failures are
terminal and prepared artifacts are checksum-recorded before alignment.
PBCCS is opt-in through `ccs_mode=pbccs` and is restricted to manifest rows
marked `read_representation=raw_subreads`; ordinary BAM and canonical CCS rows
are not converted. Normalized BED12 is the only cross-backend model contract.
`reference_metadata` may supply a JSON record containing `reference_fasta`,
`fasta`, or `path` when `reference_fasta` is not supplied. ORF products can be
published independently with `orf_mode=enabled`; DIAMOND QC consumes validated
BED12 models and publishes backend/scope-separated reports when enabled.
TransGram is an opt-in run-level backend selected with
`--model_backends transgram` (or through `all`). It consumes each validated
sorted BAM, preserves both TransGram GTF products before combining them for the
common model contract, and supports optional annotation-guided reconstruction
with `--transgram_annotation`. Its default runtime is the official project SIF
dated `20250912`; sites may override `tool_containers.transgram` with an
approved immutable local or OCI image reference.

Backend selection accepts either a JSON array in a params file or a direct
comma-separated CLI value, for example `--model_backends tama,stringtie3` or
`--model_backends all`. `backend_scopes` likewise accepts `run`, `cohort`, or
`both`. Approved inventory manifests can be processed directly; provide
`fastq_cache_dir` to resolve `classification/MD5/basename` cache entries before
falling back to acquisition.

The default Slurm fan-out is four concurrent shard tasks. Resource-sensitive
failures (exit status 137/OOM or 143/termination) retry twice, as configured in
`nextflow.config`, with linearly increased memory and wall time. Other
tool, input, and validation failures terminate immediately. TAMA additionally
uses the resource-class tiers produced by the shard planner (128/256/512 GB)
and adds 64 GB and 12 hours on each TAMA retry, up to `tama_max_retries`.
Exhausted required tasks terminate the run rather than being silently ignored.

## Legacy capability disposition

| Legacy capability | New-pipeline disposition |
|---|---|
| Manifest inspection, approved-manifest validation, path resolution | Ported into `PREPARE_LONG_READS` with typed schema/runtime validation. |
| FASTQ/BAM validation, acquisition, BAM-to-FASTQ, PBCCS | Ported as explicit preparation stages; no implicit acquisition or conversion. |
| Minimap2 alignment, samtools QC, alignment features | Ported into `ALIGN_LONG_READS`; reference and artifact identities are recorded. |
| Workload inspection, contig sharding, coverage/completeness audits | Ported into `PREPARE_CONTIG_SHARDS` with deterministic plans and explicit empty/skipped states. |
| TAMA, StringTie2, StringTie3, tmerge | Ported as independent backend subworkflows with native status/checksum products. |
| IsoQuant, FLAIR, Bambu, TransGram | Ported as optional backend subworkflows with pinned runtime references and backend-specific products. |
| Run/cohort merges and candidate comparison | Ported into merge, normalization, validation, and comparison subworkflows; run/cohort paths are isolated. |
| Evidence graph and retrospective evidence | Ported with explicit IDs/relations and post-run Parquet/publication verification. |
| ORF and DIAMOND annotation QC | Ported as optional downstream QC over validated normalized models. |
| Legacy soft-failure/recovery shortcuts and flat publication paths | Deliberately not retained; required failures terminate and outputs are namespaced. |

The default test command is:

```bash
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_all.params.json
```

Stub execution verifies wiring only. Real container-backed smoke tests and
large-input benchmarking are separate production gates.
