# Transcriptomic data inventory and manifest preparation

## Status

Proposed upstream component. This specification separates public-data discovery,
metadata resolution, representation classification, and candidate selection from
the `pipelines/long-read` reconstruction pipeline.

## Objective

Create a reusable inventory tool for answering:

> What transcriptomic data is publicly available for this species, what does
> each artifact actually contain, which samples and technologies are useful,
> and which reviewed artifacts can be handed to a downstream processing
> pipeline?

The component must support long-read data first, while keeping its core data
model general enough for short-read and other transcriptomic assays. It must
produce a complete, reviewable inventory even when no data is suitable for
processing. It must not perform transcript reconstruction, alignment, model
merging, or biological interpretation.

The output of this component is a versioned manifest and evidence package. The
long-read adapter emitted by the component is consumed by
`pipelines/long-read` through its explicit approved-manifest interface.

## Boundary with `pipelines/long-read`

The inventory component owns:

- species/taxon and repository queries;
- ENA, NCBI/SRA, BioSamples, and future provider metadata retrieval;
- artifact expansion and sidecar discovery;
- lightweight remote probing and metadata comparison;
- declared-versus-observed technology and representation classification;
- biological sample grouping and candidate selection;
- review state, quarantine state, and approval evidence;
- generation of downstream-specific candidate and approved manifests.

The processing pipeline owns:

- acquisition after an approved manifest has been selected;
- checksum verification of acquired files;
- FASTQ/BAM validation and staging;
- PBCCS conversion when explicitly requested;
- reference-aware alignment, QC, sharding, model reconstruction, merging,
  normalization, evidence collection, and final publication.

The processing pipeline must not query public archives or infer a sample's
technology from filenames. Its `discover` and `resolve_metadata` preparation
stages should either be removed from that entry point or reduced to consuming
inventory products; they must not become a second discovery implementation.

The proposed implementation belongs in a separate component, preferably
`pipelines/transcriptomic_inventory/`, with its own entry point, configuration,
schemas, provider adapters, helpers, fixtures, and tests. It may expose a
Python CLI for cache/replay and a thin Nextflow wrapper for parallel provider
work, but `pipelines/long-read/` must remain a consumer rather than an importer
of its modules.

## Product layout

An inventory run writes a self-contained directory such as:

```text
inventory/
  run_manifest.json
  inventory.tsv
  artifacts.tsv
  observations.tsv
  classifications.tsv
  selection.tsv
  candidate_manifest.tsv
  approved_manifest.tsv
  review_template.tsv
  reports/
    inventory_summary.tsv
    quarantine.tsv
    missing_metadata.tsv
    provider_status.tsv
  cache/
    provider responses and raw metadata snapshots
```

The complete inventory is retained even when only a subset is selected. The
candidate manifest is machine-proposed and is never equivalent to approval.
The approved manifest is emitted only from explicit review decisions or an
explicit, auditable safe-auto-approval mode.

Every output records:

- inventory version and tool version;
- query parameters and provider endpoints;
- taxon identifier and resolved scientific name;
- assembly/reference identity when supplied;
- UTC creation time;
- hashes of input query files, configuration, and generated manifests;
- cache policy and whether network access was required;
- provider response status and partial-failure information.

## Core inventory record

The core record is provider-neutral. Provider-specific fields may be retained
in a namespaced JSON column or sidecar, but the normalized TSV must contain at
least:

### Identity

`inventory_id`, `taxon_id`, `scientific_name`, `provider`, `run_accession`,
`sample_accession`, `experiment_accession`, `study_accession`, `biosample_id`,
`library_id`, `sample_id`, and `cohort_id` when assigned.

### Biological context

`tissue`, `organism_part`, `cell_type`, `cell_line`, `developmental_stage`,
`life_stage`, `age`, `sex`, `strain`, `breed`, `isolate`, `sample_title`, and a
free-text description. Raw provider values and normalized controlled values
must both be retained. Unknown values must be represented explicitly, not as
empty strings.

### Assay and technology

`assay_class`, `library_strategy`, `library_source`, `library_selection`,
`technology`, `platform`, `instrument_model`, `read_layout`, `strandness`,
`read_representation`, `preprocessing_state`, and `observation_unit`.

The first release must define these long-read values:

| Field | Values or examples |
|---|---|
| `assay_class` | `transcriptomic`, `other`, `unknown` |
| `technology` | `ont`, `pacbio`, `illumina`, `other`, `unknown` |
| `read_representation` | `ont_cdna`, `ont_direct_rna`, `ccs`, `hifi`, `flnc`, `raw_subreads`, `paired_end`, `single_end`, `unknown` |
| `preprocessing_state` | `processed`, `ccs_applied`, `needs_ccs`, `ccs_not_applied`, `needs_basecalling`, `needs_segmentation`, `unknown` |
| `observation_unit` | `ont_read`, `consensus_molecule`, `subread_record`, `short_read_fragment`, `unknown` |

These fields describe the evidence and processing state; they must not be
inferred solely from a filename or a repository platform label.

### Artifact and provenance

Each physical artifact is a separate row in `artifacts.tsv`, linked to the
inventory record by `inventory_id`. Required fields are:

`artifact_id`, `inventory_id`, `source_provider`, `source_uri`, `basename`,
`format`, `artifact_role`, `size_bytes`, `md5`, `sha256` when locally
available, `index_uri`, `index_format`, `index_md5`, `index_sha256`,
`submitted_or_original`, and `availability_status`.

Artifact roles include `primary_reads`, `index`, `metadata`, `annotation`, and
`unknown`. ENA semicolon-separated URI and checksum fields must be expanded
without losing their original ordering or source values.

The inventory phase may use HTTP metadata requests, range requests, and small
header probes. It must not download complete read files by default.

## Provider and discovery adapters

The provider layer must be pluggable. The initial implementation should have:

1. An ENA read-run adapter supporting `tax_eq(taxon_id)` and `tax_tree(taxon_id)`
   queries restricted to transcriptomic libraries and requested technology.
2. An ENA file-report adapter for run-level artifact and checksum metadata.
3. A BioSamples adapter for biological context enrichment.
4. An NCBI/SRA Original Format adapter for PacBio file identity and processing
   state.
5. A local-manifest adapter for replaying a frozen candidate set without
   network access.

Provider adapters must return raw responses into the cache and normalized rows
into the common model. A provider outage must be visible in `provider_status`
and must not silently convert an unresolved record into a safe candidate.

The query must be parameterized by:

- taxon ID, with optional subordinate taxa;
- assay classes and technologies;
- study, sample, or accession allow/deny lists;
- minimum metadata requirements;
- target candidate count;
- soft download and technical-risk budgets;
- probe policy and probe byte limit;
- cache directory and refresh policy.

## Evidence resolution and classification

The classification result must retain independent columns for:

1. declared provider metadata;
2. submitted artifact metadata;
3. NCBI Original Format metadata;
4. observed header/probe evidence;
5. the normalized classification;
6. confidence and reason codes;
7. the proposed downstream action.

Evidence precedence for PacBio is:

1. NCBI Original Format metadata;
2. explicit submitted artifact format and sidecars;
3. observed FASTQ header structure;
4. filename hints as weak evidence only.

The classifier must distinguish at least:

- ONT cDNA FASTQ;
- ONT direct-RNA FASTQ;
- PacBio CCS/HiFi FASTQ or BAM;
- PacBio processed/FLNC products;
- PacBio subreads BAM requiring CCS;
- PacBio raw subread FASTQ with groupable headers;
- PacBio raw subread FASTQ without reliable molecule grouping;
- mixed, conflicting, unknown, and probe-failed inputs.

FASTQ probing must record the number of records sampled, header signatures,
subread/CCS/ONT identifiers observed, malformed records, and probe errors. A
probe is evidence, not a replacement for provider metadata.

The following states are mandatory:

| State | Meaning |
|---|---|
| `READY_FOR_REVIEW` | Evidence is internally compatible and a reviewer may approve it. |
| `QUARANTINED` | Conflict, unknown representation, malformed metadata, or unsafe artifact. |
| `MIXED` | More than one incompatible representation is present for the record. |
| `CONFLICT` | Declared, submitted, original, and observed evidence disagree. |
| `PROBE_FAILED` | Classification could not be completed because the probe failed. |
| `NOT_SELECTED` | Safe enough to retain in the inventory but not chosen by the selection policy. |
| `APPROVED` | Explicit review decision exists and the downstream adapter has validated it. |

Raw PacBio subread policy must preserve these rules:

- raw subreads are not consensus molecules;
- subread base/read counts are not counted as independent transcript coverage;
- a subreads BAM may be approved for PBCCS only with the required sidecar;
- raw-subread FASTQ may be approved only through an explicit policy decision,
  with `raw_subreads` and `ccs_not_applied` retained;
- Kinnex/MAS-seq segmentation requirements remain visible;
- basecalling-required or otherwise incompatible rows remain quarantined.

## Candidate selection

Selection is a transparent ranking step, not a claim that the selected panel
is biologically complete.

The default policy must:

1. retain every discovered run in the inventory;
2. group runs by BioSample/sample identity;
3. choose at most one technical representative per biological sample for the
   first-pass panel;
4. prioritize marginal biological diversity across tissue, cell type,
   developmental stage, sex, strain, and mixed-tissue breadth;
5. use read/base counts and estimated artifact size as evidence and cost
   signals, not as proof of transcriptome quality;
6. apply technical representation, download size, and raw-subread limits as
   soft penalties by default;
7. record `SELECTED`, `SELECTED_WITH_WARNING`, or `NOT_SELECTED` plus a
   machine-readable reason;
8. allow an explicit policy to require hard exclusion for a technology or
   representation.

The selector must be deterministic for a fixed inventory and configuration.
Ties must be resolved by stable accession ordering. Its output must include
the marginal gain and every applied penalty so a reviewer can reproduce the
decision.

## Review and approval

The review file must be editable without modifying the raw inventory. It must
contain one row per candidate, the proposed action, reviewer decision, reviewer
identity, timestamp, and optional notes.

Approval rules:

- no row is approved by default;
- quarantined, conflict, mixed, or probe-failed rows require an explicit
  override reason and reviewer identity;
- approval must reference the classification report hash and selected artifact
  identity;
- a missing or empty approved set is a terminal error for downstream
  processing;
- changing a source URI, checksum, representation, or processing state
  invalidates the approval and requires re-review.

An `--auto-approve-safe` mode may approve only `READY_FOR_REVIEW` rows whose
artifact, representation, sidecar, and downstream adapter checks pass. It must
write an audit file and must never approve quarantined or ambiguous rows.

## Long-read processing adapter

The first downstream adapter emits the manifest consumed by
`pipelines/long-read` in `approved_manifest` mode. It must map inventory fields
to these processing fields:

```text
id
cohort_id
run_id
technology
assembly_id
read_set_id
sample_id
library_id
read_representation
preprocessing_state
observation_unit
source_uri
source_index_uri
input                 # present for local/replay mode; may be absent for acquire mode
minimap2_preset
minimap2_args
isoquant_data_type
expected_header_representation
classification
classification_report_sha256
```

The adapter must validate that:

- every row has one primary artifact and any required index;
- the representation and preprocessing state are compatible;
- direct-RNA rows carry their explicit alignment settings;
- raw-subread rows carry the required CCS decision;
- IsoQuant/Bambu eligibility is explicit rather than inferred from filenames;
- the selected reference/assembly identity is recorded;
- all source and classification hashes are retained.

The normal handoff is conceptually:

```bash
inventory-transcriptomic-data \
  --taxon-id 9913 \
  --assay long_read \
  --tree \
  --output-dir inventory/bos_taurus

review inventory/bos_taurus/review_template.tsv \
  > inventory/bos_taurus/approved_manifest.tsv

nextflow run pipelines/long-read/main.nf \
  --mode approved_manifest \
  --approved_manifest inventory/bos_taurus/approved_manifest.tsv \
  --preparation_stages acquire \
  --reference_fasta genome.fa
```

The exact executable name and packaging are implementation decisions; the
manifest and evidence contracts are not.

## General transcriptomic extension

Long-read classification must be an adapter over a general inventory model.
The core must also support short-read transcriptomic records with fields such
as paired/single-end layout, read length, strandedness, library preparation,
instrument model, and archive file pairs.

Future adapters may emit:

- a long-read reconstruction manifest;
- a short-read RNA-seq manifest;
- an annotation/assembly evidence manifest;
- a mixed evidence catalog for experiment planning.

The core inventory must not assume that all transcriptomic data is suitable for
long-read alignment or model reconstruction. Each adapter owns its eligibility
rules and must preserve the original inventory ID and evidence hashes.

## Reusable lessons from the legacy implementation

The following behavior is intentionally carried forward from the legacy
discovery/classification work:

- ENA taxonomic and transcriptomic filtering;
- optional subordinate-taxon searches;
- BioSamples biological-context enrichment;
- lossless provider-response caching;
- semicolon-separated artifact expansion;
- NCBI Original Format precedence over archive-renamed FASTQ names;
- lightweight FASTQ header probing;
- separation of declared platform from observed representation;
- explicit conflict, unknown, mixed, and quarantine states;
- one representative per BioSample before diversity selection;
- biology-first selection with soft download and technical penalties;
- raw-subread budgets that do not treat subread records as molecules;
- explicit sidecar, checksum, and selected-artifact provenance;
- reviewable manifests rather than implicit production approval.

The following legacy behavior is deliberately not carried forward:

- public archive discovery hidden inside the reconstruction workflow;
- a processing run that changes its own candidate set during execution;
- flat or ambiguous publication paths;
- best-effort continuation after required validation failures;
- filename-only representation inference;
- implicit conversion of raw subreads or CCS inputs.

## Proposed implementation phases

### Phase 1: extracted long-read inventory

- Move/reimplement the legacy ENA, BioSamples, NCBI Original Format, and FASTQ
  probe logic behind a standalone CLI.
- Emit the complete inventory, artifact table, classification table, selection
  table, and long-read candidate manifest.
- Add offline replay fixtures for CCS, ONT, direct-RNA, raw subreads, mixed,
  conflict, and probe-failed cases.

### Phase 2: approval and processing handoff

- Add the review-template-to-approved-manifest adapter.
- Validate the adapter against the `pipelines/long-read` manifest contract.
- Remove or disable placeholder discovery and metadata-resolution work from the
  processing entry point.
- Add an end-to-end fixture from inventory output to processing stub run.

### Phase 3: general transcriptomic inventory

- Add a provider-neutral assay model and short-read adapter.
- Add configurable provider/query plugins without changing the long-read
  manifest contract.
- Add cross-assay reports for species coverage, biological-context coverage,
  artifact availability, and unresolved metadata.

## Acceptance criteria

The inventory component is ready to replace the current placeholders when:

1. A taxon-only long-read query creates a complete cached inventory without
   downloading complete read files.
2. Replaying the same cache and configuration produces byte-stable normalized
   tables and manifests.
3. NCBI Original Format, submitted artifact, and FASTQ header evidence are
   independently visible and precedence is tested.
4. Conflicting and unknown records are quarantined and cannot enter the
   approved manifest without an override audit.
5. Selection is deterministic, sample-aware, biology-first, and fully
   explainable from its output fields.
6. The generated approved long-read manifest passes the processing pipeline's
   schema and preparation validation.
7. A processing stub run consumes the generated manifest without any import or
   runtime dependency on the inventory implementation.
8. A missing provider, missing sidecar, checksum mismatch, malformed probe,
   empty approval, and changed source artifact each fail with an explicit
   diagnostic.
9. The core inventory model can represent at least one short-read
   transcriptomic fixture without long-read-specific fields becoming required.
10. The inventory, adapter, and processing stages each publish their own
    provenance and software versions.

## Non-goals

- downloading full public datasets as part of inventory by default;
- choosing a final gene annotation or claiming transcriptome completeness;
- running minimap2, TAMA, StringTie, tmerge, IsoQuant, FLAIR, Bambu, or DIAMOND;
- silently converting input representations;
- treating a ranked candidate panel as an approved production input;
- requiring the general inventory tool to know the internals of any downstream
  reconstruction backend.
