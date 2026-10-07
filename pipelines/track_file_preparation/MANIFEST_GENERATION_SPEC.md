# Track-input manifest generation specification

## Purpose

This standalone generator runs before
`pipelines/track_file_preparation/main.nf`. It discovers input artefacts from
known filesystem layouts, resolves them into stable `run`, `sample`, or
`merged` entities, and writes the authoritative TSV consumed by Nextflow.

The generator owns discovery and entity definition. Nextflow owns track
production and validation. Nextflow must not scan production directories,
infer filename relationships, or query arbitrary input roots.

Different layouts are supported through named adapters with explicit rules,
not an expanding set of untested filename heuristics.

## Development principles

- Treat the manifest as the contract between discovery and processing.
- Resolve every path through an adapter rule, explicit pattern, or mapping file.
- Fail closed on missing required files, duplicate matches, ambiguity, assembly
  conflicts, and duplicate IDs.
- Allow `run` IDs to be derived from validated accessions, but require explicit
  grouping for `sample` and `merged` entities.
- Record source paths, matching rules, GCA, release, and configuration for
  reproducibility.
- Never modify source directories.
- Support dry-run, deterministic ordering, audit output, and exception output.
- Do not require production database access. Every output row must already have
  a GCA accession; optional UUID enrichment belongs to a later handover phase.

## Canonical output TSV

The generator writes one row per entity with this fixed column order:

```text
entity_id
entity_type
gca_accession
assembly_release
bam
gtf
track_types
bai
sj_out_tab
source_entity_ids
run_accession
sample_id
```

Production paths are absolute. Empty optional values are empty fields, not
`null`, `NA`, or fabricated paths. Tabs are forbidden in paths; commas are
reserved for `track_types` and `source_entity_ids`.

Rules:

- `entity_type` is exactly `run`, `sample`, or `merged`.
- `gca_accession` matches `GCA_[0-9]+\\.[0-9]+`.
- `track_types` is `all` or a comma-delimited subset of `coverage`,
  `gene_model`, and `splice_junction`.
- `bai` may be empty because the Nextflow pipeline creates a missing index.
- `source_entity_ids` is required for explicit sample/merged groupings.

## Source specification

The generator accepts a source-specification TSV so multiple layouts can be
combined without changing the CLI.

Required columns:

```text
source_id  source_type  root  gca_accession  assembly_release  entity_type  track_types
```

Optional columns:

```text
alignment_root  annotation_root  entity_map  bam_pattern  bai_pattern
sj_pattern  gtf_pattern  id_regex  id_replacement  sample_group  merged_id
```

Recommended CLI:

```text
--source_spec <sources.tsv>
--output_manifest <entity_manifest.tsv>
--output_audit <manifest_audit.tsv>
--output_exceptions <manifest_exceptions.tsv>
--output_summary <manifest_summary.json>
--dry_run
--strict
--default_track_types all
```

Strict mode is the production default. Exploratory mode may emit candidates,
but unresolved records must remain visible in exceptions.

## Input adapters

### `generic`

Uses explicit patterns or an entity map. Every required artefact must resolve
to exactly one path per entity. This is the fallback adapter.

### `star_run`

Supports standard STAR output:

```text
<run_dir>/alignment/<run>_Aligned.sortedByCoord.out.bam
<run_dir>/alignment/<run>_Aligned.sortedByCoord.out.bam.bai
<run_dir>/alignment/<run>_SJ.out.tab
```

It also supports the zebrafish split-root layout:

```text
<assembly>/rnaseq/output/<run>_Aligned.sortedByCoord.out.bam
<assembly>/rnaseq/output/<run>_SJ.out.tab
<assembly>/rnaseq/output/scallop/<run>_Aligned.sortedByCoord.out.gtf
```

The run ID is the BAM basename after removing
`_Aligned.sortedByCoord.out.bam`. GTF and SJ matches must use that same ID.
STAR logs, genome directories, and pass-one directories are ignored.

### `gdm_alignment`

Supports the DeepMind/GDM layout:

```text
<gdm_root>/<species>/<taxon_id>/<run_accession>/alignment/
  <run_accession>_Aligned.sortedByCoord.out.bam
  <run_accession>_Aligned.sortedByCoord.out.bam.bai
  <run_accession>_SJ.out.tab
```

The directory name and BAM basename must agree. Annotations may come from a
separate configured `annotation_root`; they must not be assumed to exist in
the alignment directory. An unreadable LTS root is an explicit access failure.

### `stringtie`

Supports files such as:

```text
<stringtie_root>/ERR2402968_2.fastq.gz.stringtie.gtf
<stringtie_root>/ERR2402972_2.fastq.gz.stringtie.gtf
<stringtie_root>/annotation.gtf
```

Filename-to-run conversion requires an explicit entity map, a configured and
tested regex/replacement, or an exact accession filename. `annotation.gtf`
becomes a `merged` entity only when explicitly declared with `merged_id`.

### `scallop`

Resolves per-run GTFs using a configured suffix, normally
`<run>_Aligned.sortedByCoord.out.gtf`, and joins them to IDs established by an
alignment adapter. Orphan GTFs are errors unless annotation-only entities are
explicitly enabled.

### `explicit_group`

Creates sample or merged rows from an explicit grouping table. Membership must
never be inferred from directory names alone.

## Resolution algorithm

For each source specification:

1. Validate source configuration, GCA, release, and adapter name.
2. Enumerate primary candidates in lexical order.
3. Derive or read the entity ID.
4. Resolve required companions using the adapter's exact rule.
5. Resolve an explicit or adjacent `.bai`/`.csi`; otherwise leave `bai` empty.
6. Check existence and readability of required paths.
7. Check entity and sanitized filename uniqueness within each GCA.
8. Emit a manifest row and audit record only after all checks pass.

No unresolved row enters the production manifest.

## Track dependencies

- `coverage` requires BAM; an absent index is valid because Nextflow creates it.
- `gene_model` requires GTF.
- `splice_junction` requires SJ unless BAM derivation is explicitly enabled
  in the downstream pipeline.
- Derivation is not silently enabled by the generator.

## Assembly and metadata

GCA precedence is: explicit source-spec value, validated path component such as
`GCA_052040795.2`, then an explicit assembly mapping table. Conflicts fail the
affected source.

Optional metadata lookup may query the assembly table by GCA and may return an
annotation UUID only when one exists. Missing UUIDs remain empty. Credentials
must not appear in manifests, audits, or command logs.

## Audit outputs

Emit alongside the manifest:

```text
manifest_audit.tsv
manifest_exceptions.tsv
manifest_summary.json
```

Audit fields should include `source_id`, `source_type`, `entity_id`,
`entity_type`, `gca_accession`, source path, resolved BAM/BAI/GTF/SJ paths,
matching rule, status, and warning.

Exception classes should include:

```text
unreadable_root
missing_bam
missing_gtf
missing_sj
ambiguous_match
duplicate_entity
assembly_mismatch
invalid_identifier
unsupported_layout
mapping_required
```

## Concrete source examples

### Zebrafish STAR + Scallop

```text
source_id= zebrafish_rnaseq
source_type= star_run
alignment_root= /hps/nobackup/flicek/ensembl/genebuild/jackt/main/zebrafish/danio_rerio/GCA_052040795.2/rnaseq/output
annotation_root= /hps/nobackup/flicek/ensembl/genebuild/jackt/main/zebrafish/danio_rerio/GCA_052040795.2/rnaseq/output/scallop
gca_accession= GCA_052040795.2
assembly_release= GCA_052040795.2
entity_type= run
track_types= all
```

This resolves the 236 currently matched BAM/GTF/SJ run entities. Missing BAI
values remain empty.

### DeepMind/GDM mouse

```text
source_id= mouse_gdm
source_type= gdm_alignment
root= /lts/production/fergal/genebuild/gdm_data/mus_musculus_domesticus_gca921998345v2/10092
gca_accession= GCA_921998345.2
assembly_release= GCA_921998345.2
entity_type= run
track_types= coverage,splice_junction
annotation_root= /path/to/approved/mouse/annotations
```

The generator must run on a host with permission to read `/lts`.

### StringTie wheat

```text
source_id= wheat_stringtie
source_type= stringtie
root= /gpfs/production/flicek/ensembl/genebuild/vianey/large_plants/25_10_aegis_wheat/GCA_965645305.1/stringtie_output
gca_accession= GCA_965645305.1
assembly_release= GCA_965645305.1
entity_type= run
track_types= gene_model
entity_map= /path/to/wheat_stringtie_entity_map.tsv
```

The mapping is required unless the filenames can be joined by an explicit
validated regex.

## Validation and test matrix

Before production output, verify readable roots, existing required files,
one-to-one companion matches, unique IDs, consistent GCA/release values,
deterministic ordering, and matching audit/manifest counts.

Fixtures must cover STAR, GDM, StringTie, Scallop, explicit grouping, separate
alignment/annotation roots, missing BAI, missing SJ, missing GTF, duplicate IDs,
ambiguous matches, GCA mismatch, unreadable roots, and mixed input types.

Recommended command shape:

```bash
python -m pytest
python -m manifest_generator \
  --source_spec test/zebrafish.sources.tsv \
  --output_manifest /tmp/entity_manifest.tsv \
  --output_audit /tmp/manifest_audit.tsv \
  --output_exceptions /tmp/manifest_exceptions.tsv \
  --output_summary /tmp/manifest_summary.json \
  --strict
```

Two runs with identical inputs must produce identical manifest rows.

## Acceptance criteria

The generator is ready when it supports named STAR, GDM, StringTie, Scallop,
generic, and explicit-group adapters; emits the canonical TSV directly
consumable by Nextflow; handles run/sample/merged semantics without guessing;
supports separate roots and missing BAI; reports missing/ambiguous inputs;
preserves absolute paths and matching-rule provenance; supports dry-run and
deterministic output; handles unreadable LTS/protected roots explicitly; and
passes the adapter/failure test matrix.
