# Track-file preparation pipeline specification

## 1. Purpose

The track-file preparation pipeline converts per-entity alignment and
annotation inputs into browser-ready track artefacts for later metadata joins.
The pipeline is driven by an externally produced TSV manifest. A separate
manifest-generation session or script owns entity definition and discovery;
this pipeline must not infer entities by scanning directories or guessing file
relationships.

The pipeline supports three entity levels:

- `run`: one input run/accession;
- `sample`: a sample-level entity assembled by the manifest producer; and
- `merged`: a merged/cohort entity assembled by the manifest producer.

For every entity, the manifest must provide a GCA assembly accession. The
production metadata database may be queried later for handover-packet
enrichment, but database access is not part of the required track-generation
path.

The initial implementation should live as a new DSL2 pipeline beneath
`pipelines/track_file_preparation/` and follow the repository's existing
Nextflow 26.04.6+, strict-shell, Singularity-first conventions.

## 2. Scope

### In scope

- Parse and validate the external entity TSV.
- Validate assembly reference inputs and GCA identity.
- Validate BAMs and generate coverage BigWigs.
- Validate GTFs and generate gene-model BigBeds.
- Validate supplied STAR junction tables.
- Optionally derive junctions from BAMs when explicitly enabled.
- Generate splice-junction BigBeds.
- Apply per-entity track selection.
- Retry appropriate infrastructure-sensitive failures and quarantine failed
  entities without stopping unrelated entities.
- Emit checksums, tool versions, provenance, exceptions, and a final track
  manifest.

### Out of scope

- Creating the input entity manifest.
- Discovering files from production directories.
- Constructing sample or merged BAM/GTF inputs.
- Assigning or discovering annotation UUIDs.
- Required production database access during track generation.
- Metadata joins or final handover-packet assembly.
- Silent coordinate or contig-name conversion between assemblies.

## 3. Input manifest contract

The canonical input is a tab-separated file with one row per entity. The
manifest producer is responsible for resolving paths and entity relationships.
The pipeline must treat the manifest as authoritative.

Required columns:

```text
entity_id
entity_type
gca_accession
assembly_release
bam
gtf
track_types
```

Optional columns:

```text
bai
sj_out_tab
source_entity_ids
run_accession
sample_id
```

### Field semantics

| Field | Semantics |
|---|---|
| `entity_id` | Stable output and metadata-join key. Must be unique within a GCA accession. |
| `entity_type` | Exactly one of `run`, `sample`, or `merged`. |
| `gca_accession` | Assembly accession, normally matching `GCA_#########.#`. Required for every row. |
| `assembly_release` | Human-readable assembly/release identifier retained in provenance. |
| `bam` | Coordinate-sorted BAM path. Required when any selected track requires BAM input. |
| `bai` | BAM index path. If absent, the pipeline may require an adjacent `<bam>.bai` or `<bam>.csi` according to explicit configuration; it must not silently create an index unless enabled. |
| `gtf` | Per-entity GTF path. Required when `gene_model` is selected. |
| `track_types` | Comma-delimited selection from `coverage`, `gene_model`, and `splice_junction`. Empty or `all` means the default selection: attempt all three. |
| `sj_out_tab` | Optional matching STAR `SJ.out.tab`. Used only for the same entity and assembly. |
| `source_entity_ids` | Comma-delimited source entities for `sample` or `merged` records. Informational and preserved in provenance. |
| `run_accession` | Optional run/accession identifier for handover and metadata joins. |
| `sample_id` | Optional sample identifier for handover and metadata joins. |

Paths must not contain tabs. Commas are reserved as delimiters for
`track_types` and `source_entity_ids`; the manifest producer must reject or
escape comma-containing paths before writing the TSV. Relative paths may be
accepted only if resolved against an explicit `--manifest_base_dir`; absolute
paths are preferred for production runs.

The pipeline must reject:

- missing required columns;
- duplicate `entity_id` plus `gca_accession` pairs;
- unsupported entity types or track names;
- missing GCA accessions;
- malformed comma-delimited values;
- missing files required by the selected tracks;
- inconsistent BAM, GTF, and SJ entity naming where a naming rule is supplied;
- assembly metadata or contig mismatches.

## 4. Reference inputs

The pipeline requires assembly reference inputs independent of the entity
manifest:

```text
chrom.sizes
contig_names/order metadata
assembly release or identifier
```

The reference contract must include a deterministic contig order. The pipeline
must validate that every input BAM, GTF, and SJ table uses legal contig names
for the selected assembly. It must not rename contigs implicitly.

The GCA accession from each manifest row must be carried through every output
record. If a future production database lookup is enabled, connection details
(`host`, `port`, `user`, and credentials supplied through the approved runtime
mechanism) belong in configuration, never in the manifest or command line
logs. Database lookup remains optional and must not be required for track
generation.

## 5. Track selection and defaults

The pipeline-level parameter `params.track_types` may override the per-row
selection. Supported values are:

```text
coverage
gene_model
splice_junction
all
```

The default is `all`, meaning the pipeline attempts all three track types for
each entity. A row-level `track_types` value may narrow the selection. The
implementation must document precedence; the recommended rule is:

1. an explicit pipeline parameter overrides the row value;
2. otherwise use the row value;
3. an empty row value means `all`.

The track types must be independently executable. Disabling `gene_model`, for
example, must not require a GTF or run any GTF process for that entity.

Junction derivation from BAM is supported but disabled by default:

```text
--derive_junctions true
```

When disabled, an entity requesting `splice_junction` must have a matching
`sj_out_tab`; otherwise it is recorded as a failed track for that entity.

## 6. Output layout and naming

The output root is `params.outdir`. Published artefacts use the following
layout:

```text
<outdir>/
  <gca_accession>/
    <entity_id>/
      coverage/
        <entity_id>.bw
      gene_model/
        <entity_id>.bb
      splice_junction/
        <entity_id>.bb
      provenance/
        <entity_id>.provenance.tsv
      exceptions/
        <entity_id>.exceptions.tsv
    track_manifest.tsv
    entity_status.tsv
    exceptions.tsv
  pipeline_info/
    execution_report.html
    execution_timeline.html
    execution_trace.txt
```

The implementation must sanitize entity IDs for filenames while retaining the
original ID in metadata. Sanitization collisions must fail preflight. The
GCA directory prevents collisions between identical entity IDs from different
assemblies.

The final `track_manifest.tsv` should contain one row per produced track, with
at least:

```text
gca_accession
assembly_release
entity_id
entity_type
track_type
track_path
source_path
status
sha256
tool_versions
normalization_parameters
```

`entity_status.tsv` must distinguish `complete`, `partial`, `failed`, and
`skipped` entities. A partial entity must never be represented as complete.

## 7. Nextflow architecture

The entry point should remain a composition boundary. The proposed structure
is:

```text
pipelines/track_file_preparation/
  main.nf
  nextflow.config
  nextflow_schema.json
  README.md
  modules/
  subworkflows/
  test/
```

The public entity channel should preserve a metadata map through all stages:

```text
tuple val(meta), path(bam), path(gtf), path(sj_out_tab)
```

The optional SJ value must have an explicit representation rather than relying
on an unstated filename convention. If separate channels are used internally,
they must be rejoined by `meta.id` and `meta.gca_accession`, never by channel
position.

The entry workflow should call `validateParameters()` and perform small
runtime checks that cannot be expressed in the schema. It should then invoke
named subworkflows approximately as follows:

```text
READ_TRACK_MANIFEST
  → PREFLIGHT_TRACK_INPUTS
  → PREPARE_COVERAGE_TRACKS
  → PREPARE_GENE_MODEL_TRACKS
  → PREPARE_SPLICE_TRACKS
  → ASSEMBLE_TRACK_SETS
  → WRITE_TRACK_PROVENANCE
  → WRITE_TRACK_MANIFEST
```

Each subworkflow should have explicit `take`, `main`, and `emit` contracts.
The entry point must not contain tool commands or long process-level logic.

## 8. Module design principles

Every module should wrap one primary tool responsibility. Independent
conversion, sorting, validation, and publication steps must remain separate.
The expected module responsibilities are:

### Input and preflight modules

- `READ_TRACK_MANIFEST`: parse the TSV into structured records.
- `VALIDATE_MANIFEST_RECORD`: validate row fields and track selection.
- `VALIDATE_REFERENCE`: validate `chrom.sizes`, assembly metadata, and contig
  order.
- `VALIDATE_ENTITY_INPUTS`: check required paths and entity naming.
- `NORMALIZE_ENTITY_METADATA`: create canonical `meta` without mutating the
  original record.

### Coverage modules

- `SAMTOOLS_QUICKCHECK`: BAM readability and index checks.
- `VALIDATE_BAM_HEADER`: header, coordinate-sort, and contig contract.
- `BAMCOVERAGE`: BAM to BigWig using configurable normalization.
- `VALIDATE_BIGWIG`: readability, non-empty content, contig, and coordinate
  checks.

### Gene-model modules

- `VALIDATE_GTF`: structured 9-column and relationship validation.
- `GTF_TO_GENEPRED`: GTF to genePred conversion.
- `GENEPRED_TO_BED12`: genePred to BED12 conversion.
- `SORT_GENE_MODELS`: deterministic BED12 sorting and normalization.
- `BED12_TO_BIGBED`: BigBed creation using `chrom.sizes`.
- `VALIDATE_GENE_BIGBED`: BigBed readability and model-field validation.

### Splice-junction modules

- `VALIDATE_STAR_JUNCTIONS`: validate supplied `SJ.out.tab`.
- `EXTRACT_JUNCTIONS`: optional BAM-derived junction extraction, preferably
  with `regtools junctions extract`.
- `STAR_JUNCTIONS_TO_BED`: apply the explicit coordinate conversion
  `start = STAR start - 1`, `end = STAR end` while preserving score and strand.
- `SORT_JUNCTION_BED`: deterministic junction BED sorting.
- `JUNCTION_BED_TO_BIGBED`: BigBed creation with custom autoSql.
- `VALIDATE_SPLICE_BIGBED`: validate intervals, scores, strand, and contigs.

### Assembly and reporting modules

- `JOIN_TRACK_RESULTS`: join by stable metadata keys and verify cardinality.
- `CHECKSUM_TRACKS`: calculate SHA-256 checksums.
- `WRITE_ENTITY_PROVENANCE`: write one structured provenance record per entity.
- `WRITE_TRACK_MANIFEST`: write the final metadata-join manifest.
- `WRITE_ENTITY_STATUS`: summarize complete, partial, failed, and skipped
  entities.

Modules must:

- use meaningful `tag` values containing GCA and entity ID;
- preserve `meta` in output tuples;
- declare all user-facing outputs explicitly;
- emit `versions.yml` where a tool is executed;
- use `task.ext.args` for configurable tool arguments;
- use version-pinned containers or environments;
- provide stubs that create every declared output;
- print diagnostics before assertions;
- run with strict shell settings inherited from the repository;
- avoid `errorStrategy 'ignore'` inside required validation modules unless the
  failure is deliberately converted into an explicit status result.

## 9. Error, retry, and quarantine policy

The pipeline is best-effort at entity level, not silent at run level.

### Retryable failures

Retry only failures likely to be transient infrastructure problems, such as:

- interrupted filesystem reads;
- temporary remote filesystem or object-store failures;
- container startup or staging failures;
- scheduler or node failures.

Retries must be bounded and configured in the pipeline configuration. Tool
syntax errors, malformed input, invalid coordinates, contig mismatches, and
failed semantic validation are not retryable.

### Entity failure handling

After the configured retry limit, an entity/track failure should be captured
with:

```text
gca_accession
entity_id
entity_type
track_type
stage
process_name
exit_status
diagnostic_path
failure_class
```

Other entities must continue. The final status must be `partial` if at least
one requested track succeeded and `failed` if no requested track succeeded.
The pipeline should complete with a non-zero exit status only for global
failures such as an invalid manifest, invalid reference, or inability to write
the output root. A policy parameter may later make incomplete entities a
release-blocking condition, but that is separate from task-level continuation.

## 10. Tool and container strategy

The preferred initial tool stack is:

| Function | Preferred tool |
|---|---|
| BAM QC | `samtools` |
| Coverage | `deepTools bamCoverage` |
| GTF validation | Repository-owned structured validator, optionally AGAT-backed |
| GTF conversion | UCSC `gtfToGenePred` |
| genePred conversion | UCSC `genePredToBed` |
| BED sorting | BEDOPS `sort-bed` or locale-stable GNU sort |
| Junction extraction | `regtools junctions extract` |
| BigWig/BigBed generation and inspection | UCSC utilities |
| Checksums | `sha256sum` |

Use Singularity/Apptainer-compatible, version-pinned BioContainers or
equivalent immutable images. Tool containers should be assigned per process in
the pipeline configuration rather than relying on host installations. The
implementation should prefer a small number of stable tool-family images when
that does not obscure versions; otherwise use one container per tool family.

Container and tool versions must be recorded in `versions.yml` and copied into
the provenance records. The pipeline must remain compatible with the root
configuration, which enables Singularity and disables Docker and Conda by
default.

## 11. Provenance requirements

For every produced track, record:

- source BAM, GTF, and/or SJ path;
- source file size and SHA-256 where feasible;
- GCA accession and assembly release;
- entity ID, entity type, and source entity IDs;
- selected track type;
- exact tool versions and container identifiers;
- command parameters, including coverage normalization;
- coordinate-conversion rules;
- reference `chrom.sizes` checksum;
- creation timestamp and pipeline revision;
- validation status and any warnings.

Provenance must be generated from structured values, not parsed from human
formatted log text.

## 12. Configuration and schema

Pipeline-specific parameters belong in
`pipelines/track_file_preparation/nextflow.config` and
`nextflow_schema.json`. Expected parameters include:

```text
--input_manifest
--manifest_base_dir
--chrom_sizes
--assembly_release
--outdir
--track_types
--derive_junctions false
--coverage_bin_size
--coverage_normalization
--max_entity_retries
--publish_failed_entities
```

The schema should express file existence, enum values, booleans, and numeric
constraints. Runtime validation should handle cross-field rules such as:

- `gene_model` requires GTF paths;
- `coverage` requires BAM and an index policy;
- `splice_junction` requires SJ input unless derivation is enabled;
- `chrom.sizes` and `assembly_release` are required together;
- an explicit pipeline `track_types` value must be supported.

## 13. Testing and verification

The implementation must include a small deterministic test profile and
representative fixture manifest. Tests should cover:

- one `run` entity;
- one `sample` entity;
- one `merged` entity;
- all tracks selected by default;
- a single disabled track type;
- supplied SJ input;
- missing SJ input with derivation disabled;
- missing SJ input with derivation enabled;
- duplicate entity IDs;
- missing GCA accession;
- invalid contig names;
- malformed GTF;
- invalid BAM sort/index state;
- one entity failing while another succeeds;
- complete, partial, failed, and skipped status generation.

Required checks for the implemented pipeline include:

```bash
nextflow lint -o concise pipelines/track_file_preparation
nextflow config pipelines/track_file_preparation/main.nf
nextflow run pipelines/track_file_preparation/main.nf \
  -stub-run -profile test \
  --input_manifest <fixture.tsv> \
  --chrom_sizes <fixture.chrom.sizes> \
  --outdir /tmp/track-file-preparation-stub
```

The implementation should inspect at least one generated `.command.sh` for
quoting, staged inputs, shell flags, and correct handling of optional files.
Stub mode proves channel wiring only; it does not replace focused validation
tests or a real-tool run using representative BAM, GTF, and reference files.

## 14. Acceptance criteria

The pipeline is ready for review when:

1. No entity discovery occurs outside the input TSV contract.
2. `run`, `sample`, and `merged` records use the same stable public channel
   contract.
3. Track types can be selected independently, with `all` as the default.
4. Junction derivation is supported but opt-in.
5. Each process has one primary responsibility and one clear output contract.
6. Global input/reference failures stop before expensive processing.
7. Entity-level failures retry only when appropriate, then produce explicit
   exception/status records while unrelated entities continue.
8. No incomplete entity is reported as complete.
9. Output paths are deterministic and include the GCA accession and entity ID.
10. Every published track has checksums, provenance, and tool versions.
11. The final track manifest is suitable for a later metadata join.
12. Stub, configuration, lint, and focused failure tests pass.

