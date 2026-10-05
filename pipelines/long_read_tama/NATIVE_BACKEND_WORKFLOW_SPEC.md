# Native long-read model-backend workflow specification

## 1. Purpose

This specification covers the four contig-sharded native merge backends:
TAMA, StringTie2, StringTie3, and tmerge. IsoQuant and FLAIR are additional
complete-BAM discovery backends and are specified separately in
`ISOQUANT_WORKFLOW_SPEC.md` and `LONG_READ_BACKEND_EVALUATION.md`; they must
not be forced through the shard-based merge topology described here.

Rework the long-read model-building pipeline so that every enabled model
backend produces a model set using its own native collapse and merge workflow.
The pipeline must make it possible to answer, for each method independently:

- what each run produced;
- what was produced from each shard;
- what was removed or merged at each native stage;
- what the final run and cohort model sets contain; and
- how those sets compare structurally and by annotation/QC evidence.

The primary use case is:

```text
--model_backend all
    -> TAMA native final model set
    -> StringTie2 native final model set
    -> StringTie3 native final model set
    -> tmerge native final model set
```

No backend's final output may be passed through another backend's merge tool
as part of the native workflow. In particular, StringTie2 and StringTie3
must not be converted to BED12 and then merged by TAMA when the purpose is to
evaluate StringTie as a method.

This is a design specification, not a request to change the scientific
thresholds, minimap2 preset, read-input classification, or TAMA parameters.

## 2. Problem statement and evidence

The current workflow has two different concepts coupled together:

1. `model_backend` selects the per-shard model builder.
2. A single global `merge_tool` selects the accession and cohort merge tool.

The current `MERGE_LONG_READ_MODELS` subworkflow receives only validated BED12
files and chooses either TAMA or tmerge for all backends. The StringTie
collapse modules do create GTF files, but `COLLAPSE_LONG_READ_MODELS` discards
those GTFs before merging. Consequently, a run with
`--model_backend all --merge_tool tama` does not produce four independent
method results: it produces four collapse candidate sets that are all forced
through the same downstream merge semantics.

This is especially important for the observed zebrafish run, where the
current common-TAMA path produced very different final counts for the same
input cohort. Those counts are useful diagnostics, but they cannot by
themselves establish which backend is better because the merge operation is
not held native to the backend being evaluated.

## 3. Definitions

- **Approved run**: one accession represented by an approved manifest row.
- **Shard**: one reference-contig BAM produced by the existing contig
  sharding stage, or the whole BAM when `shard_mode=none`.
- **Raw shard candidates**: models produced directly from one shard by a
  collapse/assembly tool, before any cross-shard merge.
- **Native run model set**: one backend's result after it has combined all
  shards belonging to one approved run using that backend's intended merge
  semantics.
- **Native cohort model set**: one backend's result after it has combined the
  native run model sets for the cohort using that backend's intended merge
  semantics.
- **Legacy common merge**: the current compatibility mode in which all
  backend candidate BEDs are merged by one selected tool.
- **Canonical representation**: BED12 for common structural validation and
  comparison. It is not necessarily the native input or native source of
  truth for every backend.
- **Provenance record**: the metadata linking a model to backend, tool
  version, stage, run, shard, input checksum, and merge parameters.

## 4. Design principles and invariants

### 4.1 Backend isolation

Every grouping key must include `backend`. No process may group, merge, or
deduplicate records from different backends in native mode.

Required grouping keys are conceptually:

```text
backend + run accession + shard
backend + run accession
backend + cohort
```

The implementation must not rely on a filename prefix alone to preserve this
invariant.

### 4.2 Preserve native evidence until native merging is complete

The workflow must preserve the native format produced by each backend until
that backend's native merge stages have finished:

- TAMA: BED12 and TAMA reports;
- StringTie2: GTF, including StringTie attributes, and derived BED12;
- StringTie3: GTF, including StringTie attributes, and derived BED12;
- tmerge: read-level GTF, tmerge output GTF/BED12, and tmerge reports.

BED12 conversion is a compatibility and comparison representation. It must
not silently replace the native input to a backend's merge operation.

### 4.3 One input, one shard plan

All backends must use the same validated alignment BAMs and the same shard
manifest for a comparison run. The native workflows may interpret those
shards differently, but the pipeline must record the identical input shard
plan and any backend-specific omissions.

### 4.4 No silent loss

Every expected shard must end in exactly one explicit state for each enabled
backend:

```text
SUCCESS
SKIPPED_BY_REQUEST
FAILED
EMPTY_OUTPUT
```

The native run and cohort reports must list missing, skipped, failed, and
successful shards. A partial result must not be presented as a complete model
set without an explicit completeness flag.

### 4.5 Determinism

Native merge input lists must be sorted by stable identifiers, not by channel
arrival order. The sorted file list and its checksum must be published beside
each merge result.

## 5. Proposed user-facing modes

Add an explicit parameter controlling merge semantics:

```text
--backend_merge_mode native|legacy_common
```

Default: `native`.

The existing `--merge_tool tama|tmerge` becomes a deprecated compatibility
parameter used only when `backend_merge_mode=legacy_common`. It must not
silently affect native mode.

Recommended validation rules:

| Configuration | Result |
|---|---|
| `model_backend=tama`, `backend_merge_mode=native` | Native TAMA workflow |
| `model_backend=stringtie2`, `native` | Native StringTie2 workflow |
| `model_backend=stringtie3`, `native` | Native StringTie3 workflow |
| `model_backend=tmerge`, `native` | Native tmerge workflow |
| `model_backend=all`, `native` | Four independent native workflows |
| any backend, `legacy_common` | Current common merge behaviour, explicitly labelled |
| `legacy_common` without `merge_tool` | Use the current default TAMA merge and emit a deprecation warning |
| native with a backend-incompatible merge option | Fail schema/runtime validation before submission |

Optional parameters:

```text
--stringtie2_merge_args
--stringtie3_merge_args
--stringtie_merge_scope run_then_cohort
--publish_native_intermediates true|false
```

The initial implementation should support only `run_then_cohort`; additional
topologies should be added only with explicit tests because merge topology can
change transcript semantics.

## 6. Native workflow matrix

### 6.1 TAMA

```text
split minimap2 BAM shards
    -> TAMA Collapse per shard
    -> TAMA Merge all successful shards for one accession
    -> TAMA Merge all accession model sets for the cohort
    -> canonical BED12 + TAMA reports
```

The existing TAMA path is the reference for the native topology. Existing
TAMA shard-retry and skip behaviour remains, but the status files must be
carried into the native run and cohort manifests.

### 6.2 StringTie2

```text
split minimap2 BAM shards
    -> StringTie2 long-read assembly per shard, native GTF
    -> StringTie2 merge of shard GTFs per accession
    -> StringTie2 merge of accession GTFs for the cohort
    -> separate derived BED12 conversion + native GTF + StringTie reports
```

The implementation must use the pinned StringTie2 container/version for both
assembly and merge. The merge inputs are GTFs, not TAMA BED12 files. The
existing `gtf_to_bed12.py` converter may be used after each native merge to
create a structurally comparable BED12, but it must not be used to feed the
StringTie merge process.

StringTie attributes such as `cov`, `TPM`, `FPKM`, and `longcov` must be
retained in the native GTF where StringTie emits them. If a merge operation
changes or removes an attribute, the change must be recorded rather than
silently reconstructed from BED12.

### 6.3 StringTie3

StringTie3 follows the same topology as StringTie2, using the pinned
StringTie3 tool and container:

```text
split minimap2 BAM shards
    -> StringTie3 long-read assembly per shard, native GTF
    -> StringTie3 merge of shard GTFs per accession
    -> StringTie3 merge of accession GTFs for the cohort
    -> separate derived BED12 conversion + native GTF + StringTie reports
```

StringTie3 must remain distinguishable from StringTie2 in process names,
version records, output paths, and reports. A common parameterized module is
acceptable only if the selected executable, image, version, and backend value
are all explicit in the task metadata and output manifest.

### 6.4 tmerge

```text
split minimap2 BAM shards
    -> BAM to read-level exon GTF per shard
    -> tmerge collapse per shard
    -> native tmerge merge of shard results per accession
    -> native tmerge merge of accession results for the cohort
    -> separate canonical BED12 conversion + tmerge reports
```

tmerge must continue to consume read-level exon GTFs. It must not consume
TAMA- or StringTie-collapsed models. The existing tmerge merge modules should
be reused or split only where necessary to make the run/cohort contracts
explicit.

## 7. Channel and module contracts

The public subworkflow contracts must carry `meta` first, as required by the
repository conventions, and must preserve backend identity explicitly. Do not
put `backend` before `meta`: that makes the contract inconsistent with the
repository's tuple handling and encourages downstream modules to unwrap
metadata differently for each backend.

### 7.1 Collapse outputs

Replace the current BED-only public contract with separate named channels or a
single documented structured tuple. Preferred shape:

```text
raw_models:
tuple val(meta), val(backend), val(shard), path(native_model), path(bed12), path(status)
```

For StringTie, `native_model` is GTF. For TAMA it is BED12. For tmerge it is
the native tmerge output representation. If a backend cannot provide one
field, emit an explicit empty/absent product through a separate named channel
rather than changing tuple arity conditionally.

Also emit:

```text
shard_audit:
tuple val(meta), val(backend), path(shard_manifest), path(status), path(stderr)
versions
```

The exact tuple syntax may differ to fit Nextflow channel constraints, but
the distinction between native files, comparison BED12, and status/audit
files must remain visible at the subworkflow boundary.

### 7.2 Native run merge inputs

Each backend-specific merge process receives only one backend and one
accession:

```text
tuple val(meta), val(backend), path(sorted_native_inputs), path(input_manifest)
```

The input manifest must contain one row per expected input, including status,
path, checksum, shard, and source process/task identifier.

### 7.3 Native merge outputs

Every backend must emit a common logical result, even though native files
differ:

```text
tuple val(meta), val(backend), path(native_model), path(comparison_bed12),
      path(model_manifest), path(model_stats), path(merge_filelist), path(versions)
```

At cohort scope, `meta` must identify the cohort and the result must include
the list of contributing accessions. A native GTF is required for StringTie2,
StringTie3, and tmerge where available; TAMA may use BED12 as its native
canonical output and derive GTF only for downstream interoperability.

### 7.4 Validation boundary

`VALIDATE_COMBINED_MODELS` must run independently for each backend's native
cohort result. It must never receive a mixed channel in which a backend label
has been lost or inferred from a filename.

## 8. Proposed subworkflow layout

Retain the existing high-level phases, but change the merge boundary:

```text
PREPARE_LONG_READS
    -> ALIGN_LONG_READS
    -> COLLAPSE_LONG_READ_MODELS
         emits raw native products and shard audits
    -> MERGE_NATIVE_LONG_READ_MODELS
         branches by backend and performs run/cohort native merges
    -> VALIDATE_COMBINED_MODELS
         once per backend
    -> RUN_DIAMOND_QC
         once per backend, with backend in metadata
```

Recommended implementation units:

- `modules/stringtie_merge.nf`: StringTie `--merge` wrapper, configurable for
  StringTie2/StringTie3 only through explicit executable/container metadata;
- `modules/native_model_manifest.nf`: deterministic file-list and checksum
  generation;
- `modules/bed12_from_native_model.nf`: conversion to comparison BED12;
- `subworkflows/merge_native_long_read_models.nf`: backend routing and
  run/cohort grouping;
- `subworkflows/merge_legacy_common_models.nf`: current common merge path,
  retained temporarily and clearly labelled;
- updates to `collapse_long_read_models.nf`: retain StringTie GTF outputs and
  expose them through the public contract;
- updates to `validate_combined_models.nf` and `run_diamond_qc.nf`: consume
  backend-labelled final outputs without special-casing TAMA as canonical.

The following boundary is mandatory, not merely an implementation preference:
conversion to comparison BED12 is a separate process from native collapse or
native merge. A StringTie process must produce GTF only; a TAMA process must
produce TAMA output only; and a tmerge process must consume and produce the
native tmerge representation. Conversion processes may run after those native
stages and may never feed a native merge process.

### 8.1 Dependency and container matrix

Every process must declare the runtime for every command in its script. The
current implementation has several failure modes this matrix is intended to
prevent:

| Process responsibility | Required commands | Required runtime rule |
|---|---|---|
| StringTie2/3 collapse | `stringtie` only | StringTie image/version matching the backend; no Python converter in the same process |
| TAMA collapse/merge | TAMA executable(s) only | One pinned TAMA image or a documented, tested split if collapse and merge need different images |
| BAM to read-level GTF | `samtools`, POSIX `awk`, `sort` | samtools image is sufficient; do not add a Python-only helper here |
| tmerge collapse/merge | `tmerge` plus native-format helpers | Image must contain tmerge and its declared Python dependencies; BED/GTF conversion is separate when it is not part of native tmerge |
| GTF↔BED12 conversion | `python3` and the staged helper | Dedicated Python image, or an explicitly verified image that contains Python; never assume a StringTie/samtools image contains it |

The `bin/` directory makes helper scripts available on `PATH`; it does not
provide their interpreter or third-party Python packages inside a container.
Each conversion module must therefore have a container smoke test that runs
`python3 --version` and the helper's `--help` (or an equivalent import test)
before it is accepted. A successful stub run does not satisfy this check.

The same rule applies to dynamically selected containers: a parameter such as
`tama_container` or `ensembl_genes_container` must have a pinned production
default, a documented architecture check, and a preflight command that
reports the image digest and executable versions. A workstation-successful
Docker image is not evidence that the image is runnable by the target
amd64 Singularity/Apptainer environment.

The entry point should remain a composition boundary. Tool commands belong in
modules; backend selection and channel grouping belong in subworkflows.

## 9. Output layout

The native output layout must prevent same-name files from different backends
or accessions overwriting one another:

```text
${outdir}/
  native_models/
    tama/
      ${accession}/
        shards/
        run_models.bed
        run_model_manifest.tsv
        run_model_stats.tsv
        cohort/
          cohort_models.bed
          cohort_model_manifest.tsv
          cohort_model_stats.tsv
          cohort_merge_filelist.tsv
    stringtie2/
      ${accession}/...
    stringtie3/
      ${accession}/...
    tmerge/
      ${accession}/...
  comparison/
    exact_and_tolerant_model_matches.tsv
    backend_summary.tsv
  qc/
    diamond/
      ${backend}/...
  reports/
    native_backend_run_report.tsv
    native_backend_cohort_report.tsv
    shard_status.tsv
```

The existing `backend_models/` path may be retained as a compatibility alias,
but new native outputs must be namespaced by backend and stage. Published
filenames must include backend/accession when a flat publish directory is
used.

## 10. Provenance and audit requirements

Every native run and cohort result must record:

- backend and native stage;
- tool name, executable version, container/image digest where available;
- pipeline revision and Nextflow version;
- reference FASTA path and checksum;
- input BAM and index checksums;
- input shard manifest and checksum;
- approved accession, sample metadata, and cohort ID;
- exact command-line arguments;
- native merge file list and checksum;
- expected, successful, skipped, failed, and empty shard counts;
- raw candidate count, post-run-merge count, and post-cohort-merge count;
- model ID mapping between native and BED12 representations;
- normalization or repair actions, including removed zero-length BED12 blocks;
- whether the result is complete or partial.

The reports must distinguish these events:

```text
read not aligned
aligned but non-primary/supplementary
aligned and present in a shard
present in a shard but not represented in a native candidate
represented in a candidate later removed/merged
present in the final model set
```

Exact read-to-final-model accounting is a separate assignment analysis and
must not be inferred from StringTie `cov`, `longcov`, or transcript counts.
The native workflow should preserve the inputs needed for that analysis and
report when it has not been performed.

## 11. Comparison and QC semantics

The pipeline must compare native final sets only after each backend has
completed its own workflow. The comparison report should include, at minimum:

- model counts at shard, run, and cohort stages;
- number of contributing accessions and successful shards;
- transcript length, exon count, intron count, and locus distributions;
- exact intron-chain matches;
- intron-chain matches with transcript-end tolerances of 25, 100, and 500 bp;
- reciprocal exon overlap and same-locus matches;
- reference/core-DB overlap when a core DB BED export is supplied;
- ORF and Diamond support using the same database and thresholds;
- per-run contribution and model retention;
- models unique to each backend under each matching rule.

Exact BED12 equality must be labelled as a strict diagnostic, not as the sole
measure of biological equivalence. Matching must be performed on normalized,
strand-aware coordinates with explicit reference contig mapping.

## 12. Core DB comparator

The classic Ensembl core model set should be treated as an external comparator
and not as another model-building backend. Add an optional, read-only export
step or companion utility that converts a core DB's `gene`, `transcript`,
`exon_transcript`, `exon`, and `seq_region` records into normalized BED12.

Required behaviour:

1. connect using an explicit user-supplied configuration, never hard-coded
   credentials;
2. query in bounded batches by stable gene/transcript IDs;
3. export transcript exon chains, strand, contig, biotype, stable IDs, and
   database identity;
4. write a checksum and export audit record;
5. never run a large unbounded production-database query from every model
   shard;
6. make the resulting BED an input to `comparison/`, not an input to native
   merging.

The core DB is a useful Ensembl-style comparator, but it must be labelled as
an independent annotation set rather than assumed to be a complete truth set.

## 13. Failure and skip policy

- A required native merge failure is terminal for that backend's affected run
  or cohort result and must produce a failure manifest. Other backends may
  continue if their inputs are independent.
- A skipped shard is allowed only when named by the backend-specific skip
  parameter. `skip_tama_shards` remains a TAMA compatibility alias; it must
  not skip StringTie or tmerge shards.
- `skip_model_shards` must use explicit `backend:accession:shard` keys, with
  an unambiguous documented parser.
- An accession with zero successful native inputs must not produce an empty
  successful run model set. It must be marked `FAILED` or `EMPTY_OUTPUT`
  according to the backend contract.
- A cohort merge must fail or be marked partial when required accession inputs
  are absent, rather than silently merging the remaining accessions.
- BED12 normalization may repair representation-only defects such as zero
  length blocks, but every repair must be counted and reported.

## 14. Configuration and schema changes

Update `nextflow.config`, `nextflow_schema.json`, README documentation, and
parameter validation together.

Required schema changes:

- add `backend_merge_mode` with enum `native|legacy_common` and default
  `native`;
- update `model_backend` description to say that `all` produces four native
  final model sets in native mode;
- mark `merge_tool` deprecated and describe its legacy-only scope;
- add StringTie merge argument parameters only if they are passed through
  safely and recorded in versions/provenance;
- add a publish-intermediates switch if storage volume requires it;
- express invalid combinations in schema where possible and enforce the
  remaining cross-field rules in a small custom validator.

Do not change the meaning of `model_backend` to mean merge tool selection.
Backend identity must remain stable throughout the workflow.

## 15. Test plan

### 15.1 Unit and contract tests

Add tests that assert:

- StringTie GTFs are retained at the collapse boundary;
- native StringTie merges receive GTFs, never BED12 from another backend;
- TAMA receives only TAMA inputs in native mode;
- tmerge receives read-level GTF/correct native inputs;
- grouping keys include backend, accession, and cohort;
- sorted merge file lists are deterministic;
- output filenames cannot collide across backend/accession;
- skipped and failed shards appear in audit outputs;
- metadata is not double-wrapped when passed to QC;
- invalid parameter combinations fail before process submission.

### 15.2 Stub matrix

Run at least:

```text
model_backend=tama,      backend_merge_mode=native
model_backend=stringtie2,backend_merge_mode=native
model_backend=stringtie3,backend_merge_mode=native
model_backend=tmerge,    backend_merge_mode=native
model_backend=all,       backend_merge_mode=native
model_backend=all,       backend_merge_mode=legacy_common, merge_tool=tama
model_backend=all,       backend_merge_mode=legacy_common, merge_tool=tmerge
```

The fixtures must contain at least two accessions and two contig shards with
deliberately different transcript structures. Assertions must verify that:

- each native backend has its own run and cohort output;
- changing `merge_tool` in legacy mode does not alter native mode;
- StringTie native output retains GTF attributes;
- the four backend outputs are not identical merely because they shared a
  merge process;
- a missing shard, failed shard, and requested skip are distinguishable.

### 15.3 Focused failure tests

Test missing/empty GTFs, malformed BED12, invalid contig names, duplicate
model IDs, incomplete accession input, checksum mismatch, non-deterministic
file ordering, and a native merge command failure. Each test must leave an
actionable status report.

## 16. Acceptance criteria

The rework is complete when all of the following are true:

1. `--model_backend all --backend_merge_mode native` produces four independent
   cohort model sets.
2. No native StringTie output is merged by TAMA or tmerge.
3. StringTie2/3 native GTFs survive through their native run and cohort merge
   stages and are published with attributes and provenance.
4. tmerge uses read-level GTF input and native tmerge merging at both scopes.
5. TAMA retains its current native collapse/merge semantics.
6. Every backend has per-shard, per-run, and cohort counts plus explicit
   skipped/failed/empty statuses.
7. Diamond and structural validation run once per backend final model set,
   with backend-safe filenames and metadata.
8. The old common merge behaviour remains available only through an explicit,
   documented compatibility mode.
9. The stub matrix and focused failure tests pass.
10. A real HPC run on the same approved manifest can be compared across
    native outputs without relying on exact BED12 equality alone.

## 17. Rollout sequence

1. Implement the new contracts and native StringTie merge modules behind
   `backend_merge_mode=native`.
2. Keep the current common merge path untouched behind
   `backend_merge_mode=legacy_common`.
3. Run the complete stub and failure matrix.
4. Run one small approved manifest with all four native backends.
5. Run the full cohort with the same reference, alignment, shard manifest,
   Diamond database, and resource profile used for the current comparison.
6. Generate the native per-run report and compare it with the legacy report.
7. Inspect StringTie2/3 retention, native merge counts, QC support, and core
   DB overlap before selecting a production default beyond native mode.
8. Deprecate, but do not immediately remove, `legacy_common` after one release
   cycle and after existing result consumers have migrated.

## 18. Preflight gates before implementation is accepted

The redesign must not be called ready because the schema parses or a stub
workflow completes. The following gates are required in order:

1. **Contract gate** — the schema, pipeline defaults, README, and runtime
   validation agree on mode names, defaults, legacy-only parameters, and
   conditional requirements. In particular, `native` must not become the
   default until the entry point actually dispatches to native workflows.
2. **Boundary gate** — every public tuple is `meta` first; backend, accession,
   cohort, stage, and native/comparison representation are explicit; no
   filename is used as backend identity.
3. **Dependency gate** — each process has one primary tool responsibility and
   every command in its script exists in its declared image. This includes
   helper interpreters and Python packages, not only the headline tool.
4. **Container gate** — every literal image is reachable and amd64-compatible;
   dynamic images have pinned defaults and an explicit inspection path. A
   failed network check remains unverified rather than being treated as a
   passing image check.
5. **Failure gate** — native merge and validation failures are terminal for
   the affected backend/stage. Global `errorStrategy 'ignore'` must not turn a
   required native failure into a successful partial result.
6. **Mode gate** — run the full native/legacy matrix, including one invalid
   combination that must fail before process submission and one missing,
   failed, empty, and requested-skip shard case.
7. **Execution gate** — lint, config parse, schema validation, stub matrix,
   focused failure tests, and at least one real container smoke test have
   passed. Stub success is wiring evidence only.

The native implementation now satisfies the basic boundary and mode gates:
`COLLAPSE_LONG_READ_MODELS` emits native-format products, backend-specific
run/cohort merge processes are selected in native mode, and GTF-to-BED12
conversion is separate. Expected/skipped/failed shard manifests are now
emitted for every requested backend, native merge arguments and file-list
checksums are recorded, and required native collapse/merge failures are
terminal. The remaining acceptance work is operational: validate the real
containers on an amd64 execution host and perform representative real-input
smoke/HPC execution before treating the redesign as production-ready.
