# Long-read implementation handoff

## Project-level objective

Replace `pipelines/long_read_tama/` with a strict, contract-driven DSL2
long-read reconstruction pipeline under `pipelines/long-read/`, covering
explicit preparation, reference-aware alignment/QC, deterministic workload and
contig sharding, selectable reconstruction backends, native merges, validated
BED12 normalization, optional evidence/ORF/DIAMOND products, and reproducible
provenance. The new tree has no runtime import from the legacy tree.

## Current hypothesis / mental model

The pipeline is composed as `PREPARE_LONG_READS → ALIGN_LONG_READS →
PREPARE_CONTIG_SHARDS → selected backend subworkflows →
MERGE_LONG_READ_MODELS → NORMALIZE_VALIDATE_MODELS → optional evidence/QC →
FINALIZE_LONG_READ_PIPELINE`. `meta` is retained on public tuples; plans and
checksums are authoritative provenance, while BED12 is the only cross-backend
model contract.

## Invariants and decision rules

* Required stages terminate on failure; there is no `errorStrategy 'ignore'`.
* `all` expands deterministically to the seven supported backends.
* Empty contigs remain in the shard plan and do not launch split tasks.
* Raw PacBio subreads are the only input representation eligible for PBCCS;
  canonical CCS and ordinary BAM inputs are never silently converted.
* Source/reference identity is checked before sharding and mapped-read totals
  are conserved across emitted shards.
* Shard plans classify eligible contigs as `small`, `large`, or `very_large`;
  `shard_max_concurrent` bounds shard/backend process selectors when set.
* Native outputs are never compared or published as common models without GTF
  to BED12 normalization and validation.
* The test profile disables containers only for stub wiring; production defaults
  keep Singularity enabled and `use_containers=true`.
* The legacy `pipelines/long_read_tama/` tree is source material only.
* Reusable backend modules read optional `task.ext.args` values and retain
  validated pipeline-parameter fallbacks; pipeline-owned flags remain built by
  the workflow and rejected from free-form parameter fields.

## Parameter/schema matrix

The executable contracts are `nextflow_schema.json` for the main entry point
and `retrospective_schema.json` for retrospective evidence collection.

| Group | Main-schema parameters |
|---|---|
| Input/preparation | `mode`, `manifest`, `approved_manifest`, `reference_fasta`, `reference_metadata`, `reference_fai`, `preparation_stages`, `discovery_query` |
| Identity/alignment | `outdir`, `tracedir`, `cohort_id`, `minimap2_preset`, `secondary_mode`, `alignment_qc`, `use_containers` |
| Sharding/publication | `shard_mode`, `shard_contig_reads`, `shard_very_large_reads`, `shard_max_concurrent`, `publish_native_intermediates` |
| Backend selection/scope | `model_backends`, `backend_scopes` |
| TAMA/StringTie/tmerge | `tama_cap_mode`, `tama_end_mode`, `tama_run_mode`, `stringtie2_merge_args`, `stringtie3_merge_args`, `tmerge_args` |
| IsoQuant | `isoquant_mode`, `isoquant_scope`, `isoquant_data_type`, `isoquant_genedb`, `isoquant_analysis`, `isoquant_high_memory`, `isoquant_args`, `isoquant_large_output` |
| FLAIR/Bambu | `flair_args`, `flair_threads`, `bambu_scope`, `bambu_track_reads`, `bambu_quantify` |
| Evidence/annotation QC | `evidence_mode`, `orf_mode`, `diamond_mode`, `diamond_database`, `diamond_reference_proteins` |
| Runtime registry | `tool_containers` |

The main schema contains typed enums, file-path checks, uniqueness constraints,
immutable SHA-256 container-reference patterns, and conditional
reference-guided/DIAMOND requirements. The retrospective schema declares the
inherited selector fields needed to validate the retrospective entry point
without scheduling main-pipeline work and applies the same immutable container
reference contract.

## Completed milestone and evidence

Implemented the new directory structure, typed schemas, preparation/alignment/
sharding/backend/merge/normalization/evidence/finalization subworkflows,
checked-in helpers, deterministic fixtures, channel-contract documentation,
container registry, and unit tests.

Passing checks:

* `pytest -q pipelines/long-read/test` — 29 passed, including indexed-BAM,
  explicit PBCCS sidecar, shard coverage, backend shard auditing, reference
  compatibility, alignment extraction, and BED12 checks.
* `python3 -m compileall -q pipelines/long-read/bin pipelines/long-read/test` — passed.
* `nextflow lint -o concise pipelines/long-read` — no lint errors (warnings are
  unused workflow parameters and helper staging advisories).
* `nextflow config pipelines/long-read/main.nf` and the corresponding
  retrospective command — passed.
* TAMA stub fixture — passed.
* FASTQ stub fixture — passed.
* `shard_mode=none` stub fixture — passed.
* evidence-enabled stub fixture — passed, including `catalog.json`,
  `relations.parquet`, and `verification.json`.
* DIAMOND-enabled stub fixture — passed.
* Retrospective evidence stub fixture — passed.
* all-backend stub fixture with `run,cohort` scopes — passed for TAMA,
  StringTie2, StringTie3, tmerge, IsoQuant, FLAIR, and Bambu.
* Two-run/all-backend stub fixture — passed with independent `per_run/run1` and
  `per_run/run2` model trees plus one cohort tree; native and finalization input
  filename collisions were fixed and retested.
* invalid approved-manifest parameter combination — rejected by nf-schema before
  tasks were submitted.
* Real host-tool tmerge fixture — passed after deterministic coordinate sorting
  of the converted read-level GTF; this exercised preparation, samtools
  alignment validation, shard coverage, tmerge, merge, and BED12 validation.
* Real zero-mapped-BAM tmerge fixture — passed with no split/tool task,
  explicit EMPTY native product, audited plan, EMPTY merge, and valid empty
  BED12 publication.
* Native-publication stub fixture — passed; native products were namespaced by
  run/backend/shard and `versions.yml` was excluded from native folders; the
  publication tree now follows `models/<backend>/normalized` and
  `models/<backend>/native/shards/<shard>`.
* Real host-tool tmerge provenance fixture — passed after adding reference
  identity, run/cohort IDs, normalized parameter configuration, preparation
  checksum, shard-plan checksums, and final-product checksums to the run
  manifest.
* All production tool references are immutable registry digests in both
  `nextflow.config` and `resources/containers.yml`.
* FLAIR combined native output is now audited and used as the run-level merge
  input, while per-shard outputs remain the completeness-audit source.
* Evidence package stub and real host-tool tmerge fixture — passed with
  explicit reference/alignment/model/assessment IDs, compatibility metadata,
  checksums, a genuine Parquet relation table, and verifier output. The
  two-run evidence stub produced two alignment artifacts, 21 model artifacts,
  and 28 explicit relations. Apache Parquet CLI `meta`, `head`, and
  `check-stats` validation was also run against the generated relation table;
  all three succeeded.
* Retrospective evidence fixture — passed with a historical BAM/BED12 manifest,
  explicit alignment/model/assessment records, verified Parquet relations, and
  finalization outputs.
* The all-backend fixture now publishes a cohort `model_comparison` table,
  complete shard status/coverage/backend-audit products, and
  `inputs/reference_identity.json`.
* The current post-refactor all-backend stub completed successfully at
  `/tmp/long-read-all-tty`, producing 14 validated BED12 artifacts, one cohort
  model-comparison table, and a final run manifest through the new audit and
  comparison subworkflows.
* Retrospective validation warnings from inherited main-pipeline defaults were
  removed by declaring those inherited fields in the retrospective schema;
  only the inherited entrypoint-specific selector and macOS runtime-metrics
  warnings remain.
* Explicit preparation stages are selectable through `preparation_stages`:
  discovery, metadata resolution, acquisition, and BAM-to-FASTQ conversion.
  The combined preparation-stage stub passed with a finalized checksum-bearing
  manifest; the direct-input regression still passes.
* `reference_metadata` resolves a relative or absolute FASTA path from a
  validated JSON record, and an ORF-only stub run publishes namespaced ORF
  products without requiring DIAMOND.
* IsoQuant selected read-level products are retained under a dedicated QC
  namespace, with content checksums for reference/BAM/GeneDB inputs in a
  backend manifest. FLAIR now preserves discovered native GTF/BED/FA/read-map
  files and their checksums in per-shard/merged manifests, including explicit
  EMPTY-shard manifests.
* Bambu input provenance now records SHA-256 checksums for the BAM and
  reference alongside explicit read-tracking and quantification flags.
* Main and retrospective entry points now compose named subworkflows only;
  native-shard auditing, normalized-model comparison, retrospective manifest
  validation, and retrospective evidence collection are encapsulated behind
  subworkflow boundaries rather than direct process imports.
* The README records explicit legacy capability port/adapt decisions, and the
  channel-contract document lists named subworkflow boundaries with their
  public inputs/outputs. Backend modules read optional `task.ext.args` values
  with validated parameter fallbacks.
* Merge file lists, file-list checksums, merge status, skipped-input status,
  and merge manifests are now published under namespaced run/cohort/backend
  `models/<backend>/merge/` directories. Merge manifests now embed the sorted
  source input list, file-list checksum, backend, scope, and status; merged GTF
  intermediates remain private unless separately selected for publication. The
  final run manifest also carries the normalized merge input lists and their
  file-list SHA-256 values.
* Finalization now runs `bin/verify_publication_checksums.py` and publishes a
  checksum-verification report; the helper has tamper-detection and published-
  path regression tests. Final checksum rows identify the eventual published
  model path separately from the staged source path used for verification. The
  current real host-tool tmerge fixture produced zero checksum errors and
  checksums for the published BED12, model manifest, model summary, and
  validation report, with no task-internal `.command.*` entries. The same
  manifest now records explicit preparation, alignment, sharding,
  reconstruction-by-backend, merge, normalization, validation, optional-QC,
  evidence, and finalization statuses; enabled evidence is a dependency of
  finalization.
* Real Docker-backed tmerge, FLAIR, IsoQuant, Bambu, and StringTie2 fixture
  runs completed on the tiny indexed-BAM fixture with immutable
  `linux/amd64` images. Each produced a final run manifest, publication
  checksums, and checksum-verification report. The optional backends returned
  explicit EMPTY products where the tiny fixture had no discoverable models.
  After correcting native-status propagation, fresh IsoQuant and Bambu Docker
  runs again completed preparation, alignment, sharding, reconstruction,
  normalization, validation, and finalization, with run- and cohort-level
  validated BED12 products under `/tmp/long-read-real-docker-isoquant` and
  `/tmp/long-read-real-docker-bambu`.
* A real Docker-backed two-run tmerge fixture completed with
  `succeededCount=36; failedCount=0`; it produced independent `run1` and
  `run2` alignment, shard, merge, normalized BED12, validation, and checksum
  products under `/tmp/long-read-real-docker-two-run/per_run/`.
* A real Docker-backed TAMA fixture using the official image's bundled
  `gmap_test.sam`/`test_genome.fa` data completed through the full pipeline,
  including SAM preparation, TAMA collapse, merge, BED12 normalization,
  validation, and finalization. The tiny repository fixture's TAMA crash is
  therefore recorded as a fixture/tool edge case, not converted into EMPTY.
* Docker container smoke passes 12/12 configured images, including the Python
  helper image, using `--engine docker --platform linux/amd64`.
* The post-refactor all-backend stub, FLAIR/IsoQuant/Bambu product stubs, and
  retrospective evidence stub were rerun after the latest process changes and
  completed successfully. StringTie3's executable/version smoke now passes in
  both Docker and Apptainer, but the real workflow fixture remains unproven on
  this Apple Silicon/QEMU host because the pinned executable previously
  segfaulted during model reconstruction under emulation. TAMA reaches the
  official tool but the tiny repository fixture triggers an upstream
  `IndexError` during locus processing.
* Invalid parameter checks were exercised with fresh stub invocations:
  reference-guided IsoQuant without `isoquant_genedb`, DIAMOND without a
  database/protein FASTA, and duplicate backend selection all failed during
  parameter validation before task submission.
* Software-version aggregation now preserves complete per-process YAML records
  instead of sorting their lines globally. A fresh real tmerge run contains 24
  process records, resolved minimap2/samtools versions, the merge backend
  record, and no unresolved command substitutions. Evidence packaging and ORF
  extraction/prediction versions are also included when those stages are
  enabled.
* The repository-local immutable-reference audit
  (`bin/audit_containers.py`) passed all twelve configured images for a
  linux/amd64 OCI descriptor or single manifest. The stale tmerge digest was
  replaced with the current Perl 5.38-slim index digest, expressed in
  digest-only form so Apptainer accepts it. After placing Apptainer scratch
  space on the Lima guest's root disk and allowing up to 1200 seconds per
  image, the checked-in smoke helper passed all 12/12 executable/version
  checks, including tmerge Perl 5.38.5, IsoQuant 3.9.0, FLAIR 3.0.1, Bambu
  3.12.1, DIAMOND 2.1.9, and the Python 3.12.16 helper. The JSON evidence is
  `/tmp/long-read-container-smoke-apptainer-1200.json`. This is complete
  Linux/amd64 smoke evidence inside the local x86_64 Lima guest, but not a
  substitute for native Linux/amd64 production execution because the guest
  itself is QEMU-emulated on Apple Silicon.
* Checked-in production-gate helpers now make the remaining external checks
  reproducible: `bin/container_smoke.py` runs the executable/version matrix
  and records JSON results, while `bin/record_benchmark.py` records the exact
  Nextflow command, wall time, trace task/retry summary, output checksums, and
  software-version artifact for a supplied historical dataset.
  The timeout/error path is covered, kills the local process group, performs
  best-effort cleanup of matching Lima guest processes, and writes a complete
  result record even when Apptainer cannot pull or start an image.
* Non-trivial pipeline-owned Python helpers used for manifest normalization,
  shard auditing, BED12 conversion/validation, evidence packaging, comparison,
  and finalization now run in the pinned `python_helper` image when containers
  are enabled. A fresh Docker tmerge fixture completed all stages with these
  helper containers and published validated models under
  `/tmp/long-read-real-docker-python-helper-tmerge`.
* `bin/record_benchmark.py` was exercised with a typed params file on the tiny
  fixture and produced a successful `benchmark.json`, trace summary, command
  record, and checksum table under
  `/tmp/long-read-record-benchmark-fixture-typed`. This validates the recorder
  mechanics only; it is not the required large historical benchmark.

## Files touched

All implementation files are new under `pipelines/long-read/`: entry points,
configuration and schemas, `modules/`, `subworkflows/`, `bin/`, `resources/`,
`scripts/`, `test/`, and `docs/`. Existing unrelated untracked files elsewhere
in the worktree were not modified.

## Remaining tasks and acceptance criteria

1. Run the complete linux/amd64 Apptainer runtime audit and record executable
   versions for every configured image on a native Linux/amd64 host. The
   Docker matrix and local extended Apptainer smoke are both 12/12, but the
   local guest is QEMU-emulated on Apple Silicon and therefore does not close
   the native production-runtime gate.
2. Re-run the real container-backed StringTie3 fixture on a native Linux/amd64
   runner or approved larger fixture. The pinned executable currently
   segfaults under Apple/QEMU emulation; do not silently substitute another
   backend.
3. Benchmark at least one large historical BAM and BED12 set and record
   resource use, task retries, output checksums, and software versions using
   `bin/record_benchmark.py`.
4. Complete the remaining backend-specific process parity review and record
   the final multi-run publication tree, including the new comparison and
   completeness reports. Docker fixtures for FLAIR, Bambu, and IsoQuant now
   pass on the tiny fixture; native Linux coverage is still required for the
   production gate.

Representative two-run/all-backend stub output tree (verified at
`/tmp/long-read-final-two-run.QYyeZM`) contains independent
`per_run/run1/`, `per_run/run2/`, and `per_cohort/fixture/` roots, each with
backend-separated normalized models. Run roots also contain alignment, QC,
transcript-feature, and shard products; the cohort root contains the reference
and model-comparison products.

The implementation is therefore wiring-complete for stub validation, but not a
production-readiness declaration until these acceptance criteria are met.

Current handoff state: repository `/Users/jackt/projects/ensembl-genes-nf`,
branch `feature/minimap_tama_hg002`, HEAD `d869d4bd`; no commit has been made
for this implementation and `pipelines/long-read/` is currently untracked.
Unrelated worktree changes were not modified.

## Non-goals and traps to avoid

Do not reintroduce legacy soft-failure behavior, infer transcript strand from
alignment orientation, match evidence by filenames, enable runtime package
installation/network access, publish every shard BAM by default, or call a stub
run biological validation. Do not treat a zero-read plan as a tool failure, but
do fail if a non-zero source BAM has no eligible or missing shard. Do not stage
`work/`, `.nextflow*`, caches, generated results, or `__pycache__` directories.

## Exact commands/tests

```bash
git diff --check
nextflow lint -o concise pipelines/long-read
nextflow config pipelines/long-read/main.nf >/tmp/long-read-main.config
nextflow config pipelines/long-read/retrospective.nf >/tmp/long-read-retrospective.config
pytest -q pipelines/long-read/test
python3 -m compileall -q pipelines/long-read/bin pipelines/long-read/test
python3 pipelines/long-read/bin/audit_containers.py \
  pipelines/long-read
python3 pipelines/long-read/bin/container_smoke.py \
  pipelines/long-read/resources/containers.yml \
  --engine docker --output /path/to/container-smoke.json --platform linux/amd64
# For Apptainer, set APPTAINER_TMPDIR to a filesystem larger than /tmp and
# allow enough time for OCI-to-SIF conversion on constrained emulated hosts:
APPTAINER_TMPDIR=/var/tmp python3 pipelines/long-read/bin/container_smoke.py \
  pipelines/long-read/resources/containers.yml \
  --engine apptainer --output /path/to/apptainer-smoke.json \
  --platform linux/amd64 --timeout-seconds 1200
# Add --disable-cache when a fresh OCI-to-SIF conversion is required.
# The generic skill helper below does not parse @sha256 references; its failure
# is diagnostic only, not the container-audit acceptance check.
python3 .codex/skills/ensembl-nf-validate-pipeline/scripts/validate_containers.py \
  pipelines/long-read --platform linux/amd64
nextflow run pipelines/long-read/main.nf -stub-run -profile test
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_all.params.json
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_publish.params.json
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_preparation.params.json
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_reference_metadata.params.json
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_orf.params.json
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_isoquant_products.params.json
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_flair_products.params.json
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_bambu.params.json

# Historical-data benchmark (run on Linux with the approved large dataset)
python3 pipelines/long-read/bin/record_benchmark.py \
  --pipeline pipelines/long-read/main.nf \
  --manifest /path/to/historical_manifest.tsv \
  --reference /path/to/reference.fa \
  --outdir /path/to/long-read-benchmark \
  --params-file /path/to/benchmark.params.json \
  --extra -ansi-log false
nextflow run pipelines/long-read/retrospective.nf -stub-run -profile test \
  --retrospective_manifest pipelines/long-read/test/data/manifest.tsv \
  --cohort_id fixture --evidence_mode enabled

# Real host-tool smoke (samtools/perl available on the developer host)
nextflow run pipelines/long-read/main.nf -profile test \
  -params-file pipelines/long-read/test/data/real_tmerge.params.json

# Multi-run publication topology
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_all.params.json \
  --manifest pipelines/long-read/test/data/manifest_two_runs.tsv \
  --outdir /tmp/long-read-two-run

# Real Docker-backed multi-run tmerge fixture
nextflow run pipelines/long-read/main.nf \
  -c /tmp/long-read-docker.config \
  -params-file /tmp/long-read-real-docker.params.json \
  --manifest pipelines/long-read/test/data/manifest_two_runs.tsv \
  --outdir /tmp/long-read-real-docker-two-run -resume -ansi-log false
```
