# TAMA workload partitioning and resource escalation specification

## Status

Proposed. The first implementation milestone is limited to handling large
whole-accession and contig-level TAMA tasks with appropriate resources.
Adaptive shard planning and coordinate-window sharding are explicitly
deferred.

## Objective

Make the long-read TAMA workflow robust for accessions whose alignment BAMs or
individual reference contigs are too large for a single default TAMA Collapse
task, while preserving the repository's canonical model-reduction order:

```text
reads for one accession
    -> TAMA Collapse
    -> optional per-accession shard merge
    -> one model set per accession
    -> final TAMA Merge across accessions
    -> one combined model set
```

The workflow must never trade memory management for silent transcript
truncation, duplicate models, or competing final model sets.

## Background and problem

TAMA Collapse reduces redundant read alignments to transcript models. It needs
the complete evidence assigned to its input scope. TAMA Merge operates later on
already-collapsed model BED files and removes/reconciles redundant models.

The current workflow has two relevant modes:

* `whole`: one TAMA Collapse over the complete sorted BAM for each accession;
* `contig`: one TAMA Collapse per reference contig, followed by one TAMA Merge
  per accession.

The `whole` mode gives simple channel topology but concentrates all memory use
in one task. A scheduler exit status of 137, such as the failure observed for
`SRR32588715:whole`, indicates that the task was killed for exceeding its
resource allocation. TAMA's `-rm low_mem` argument does not guarantee that the
process will fit within the scheduler allocation.

The safe first response is reference-aware partitioning plus resource
escalation. Arbitrary read sharding and coordinate windows are not equivalent:
they can separate evidence for one transcript or cut a transcript at a window
boundary.

## Invariants

The implementation must preserve these invariants:

1. A TAMA Collapse input contains reads from exactly one approved accession.
2. Each alignment is assigned to at most one collapse shard.
3. Each reference contig is assigned to exactly one shard when contig mode is
   used. Zero-length `idxstats` entries may be omitted.
4. A contig-sharded accession is merged back to exactly one accession model set
   before the cohort merge.
5. The final cohort merge receives exactly one model set per accession.
6. Inputs to both accession and cohort merges are deterministic and published.
7. No coordinate-window sharding is enabled by this milestone.
8. A failed validation, collapse, or merge terminates the required production
   path; retries may increase resources but must not hide tool failures.

## Phased implementation

### Milestone 1: resource escalation for oversized tasks

This milestone does not change the current sharding topology. It makes the
workload of existing whole-accession or contig-level tasks explicit, gives
large tasks larger resources, and provides a safe operational path for
accessions that cannot be processed in the current allocation.

#### Required behaviour

1. For contig-level execution, inspect the sorted BAM with its index before
   materialising shards. For whole-accession execution, retain the existing
   accession-level task.
2. Record, where available, per reference contig:
   * contig name and reference length;
   * mapped alignment count from `samtools idxstats`;
   * whether the contig is eligible for processing;
   * the selected resource class.
3. Assign each existing TAMA task to a resource class using its accession or
   contig workload metadata.
4. In contig mode, run one TAMA Collapse per eligible contig and merge all
   resulting
   BEDs once per accession.
5. In whole mode, preserve the existing one-collapse-per-accession contract.
6. Retain the current explicit `whole` and `contig` modes. Do not introduce an
   adaptive mode as part of this milestone.

#### Resource policy

Resource selection must be pipeline-local and configurable. It must not
require users to edit the shared root configuration.

The initial policy should use workload tiers rather than pretending that a
read-count target is a hard memory bound:

| Workload | Initial action |
| --- | --- |
| ordinary whole accession | normal high-memory TAMA label |
| large accession | high-memory TAMA label |
| very large contig | larger high-memory allocation |
| exit 137/140/143 | retry with the next configured allocation |
| exhausted retries | fail with accession, contig, requested memory, and task work directory |

Mapped-read count is an initial workload proxy, not a scientific threshold.
Resource tiers must be easy to override with `task.ext` or pipeline-local
configuration, and the selected estimate/class must be visible in task tags,
reports, or both.

The retry policy must be bounded. Increasing memory is appropriate for an
infrastructure kill; it must not retry malformed BAMs, invalid references, or
TAMA validation failures.

#### Parameters and nf-schema

The following are proposed public parameters. Names may be adjusted to match
existing repository conventions, but each must be represented in the
pipeline's `nextflow_schema.json` and documented in `README.md`:

* `shard_mode`: the existing enum `whole`, `contig`; any future `adaptive`
  value requires the milestone 2 design;
* `shard_contig_reads`: mapped-read threshold used to classify a contig as
  large;
* `tama_memory_small`, `tama_memory_large`, and
  `tama_memory_very_large`, or an equivalent pipeline-local resource policy;
* `tama_max_retries`: bounded resource-escalation retry count.

Defaults must be conservative and documented as operational heuristics. The
schema should validate types, non-negative values, and enum membership. Any
rule depending on BAM contents must remain a small runtime validation because
nf-schema cannot inspect a remote/indexed BAM.

The entry point must call `validateParameters()` before launching work. No
conditional process output may be accessed unless that process was invoked;
absent optional branches must use empty channels and explicit branch-safe
composition.

### Milestone 2: adaptive shard planning and grouping

Only after milestone 1 is stable, add an `adaptive` mode that selects between
whole-accession and contig-level execution. It should then add
workload-balanced grouping for small contigs. Large contigs remain
whole-contig shards. Small contigs may be packed into groups using a
deterministic greedy algorithm sorted by descending estimated workload.

The splitter must publish a shard manifest containing at least:

```text
accession  shard_id  contigs  mapped_reads  reference_bases
```

Each shard must have an indexed BAM and must retain metadata through
`TAMA_COLLAPSE`, validation, and the accession-level `TAMA_MERGE`.

The grouping target is a soft packing target. A single large contig is never
split merely to satisfy it.

### Milestone 3: coordinate-window sharding, if still necessary

Coordinate windows are out of scope for the first two milestones. They require
an explicit design for:

* window overlap width;
* preservation of reads and supplementary alignments;
* duplicate model detection across adjacent windows;
* reconstruction of models crossing boundaries;
* validation against an unsplit whole-contig result;
* deterministic provenance for models assembled from multiple windows.

No implementation may expose window sharding as a production option until
those rules and tests exist. TAMA Merge alone must not be assumed to repair
truncated boundary models.

## Proposed module boundaries

Keep orchestration in the long-read subworkflow and tool invocations in
pipeline-local modules:

* `INSPECT_BAM_WORKLOAD`: reads the BAM index and emits deterministic contig
  workload metadata;
* `PLAN_TAMA_SHARDS`: creates the selected `whole` or `contig` plan;
* `SPLIT_BAM_BY_CONTIG`: materialises indexed BAM shards for the plan;
* `TAMA_COLLAPSE`: receives one accession and one planned shard;
* `TAMA_MERGE`: merges shard BEDs per accession and then accession model sets
  for the cohort.

The public subworkflow contract should remain a tuple beginning with `meta`.
Shard metadata must include an accession identifier and stable shard ID; it
must not rely on filenames alone.

## Accuracy and completeness checks

The implementation must verify:

* BAM and BAI exist and are readable;
* every eligible contig is assigned once;
* the shard plan is deterministic for identical BAM/index inputs;
* mapped-read totals before and after extraction agree within the defined
  `samtools` accounting rules;
* every shard produces a non-empty, valid TAMA BED when reads are present;
* every accession contributes exactly one model set to the final merge;
* merge file lists are sorted and checksummed;
* no accession is silently omitted when one shard fails.

For a representative accession, compare whole and contig-sharded outputs by
structural model signature rather than model ID:

```text
reference contig, strand, transcript start, transcript end, exon coordinates
```

The comparison should report model counts, exon/junction counts, transcript
length distributions, and models present only in one topology. Identifier
differences alone are not evidence of a biological difference.

## Test plan and acceptance criteria

### Structural tests

* nf-schema accepts `whole` and `contig` and rejects invalid enum values;
* invalid thresholds and retry counts are rejected before task submission;
* `whole` and `contig` branches have no undefined output-channel access;
* the stub path produces the expected per-accession and final merge topology;
* the shard manifest and merge file lists are deterministic.

### Focused functional tests

Use a tiny indexed BAM containing multiple contigs to test:

* one contig per shard;
* omission of zero-length contigs;
* preservation of mapped-read counts;
* accession-level merge after multiple collapse outputs;
* final cohort merge from exactly one model set per accession;
* clear failure for a missing BAM index or missing shard output.

### Resource-policy tests

Use synthetic workload metadata to verify that:

* ordinary workloads select the normal tier;
* large workloads select the larger tier;
* exit 137/140/143 advances to the next tier only up to the retry limit;
* validation and biological/tool errors do not trigger resource retries.

### Representative run

Run one small real container-backed example with a tiny reference, indexed BAM,
and two accessions. Then run the same example in whole and contig modes and
compare structural signatures. A stub run or Nextflow lint result alone is not
acceptance evidence for the TAMA behaviour.

## Non-goals and traps

This work does not change TAMA end/cap parameters, minimap2 presets, read
classification, CCS generation, or model-QC thresholds. It must not:

* shard by arbitrary read batches;
* split genomic coordinates without boundary reconciliation;
* merge runs before each run has a complete model set;
* create separate competing final model sets by tissue or platform;
* use an unbounded retry loop for memory failures;
* make a read-count threshold appear to be a guaranteed memory requirement.

## Definition of done for milestone 1

Milestone 1 is complete when existing whole-accession and contig-level TAMA
tasks receive an appropriate resource tier, memory-killed tasks advance
through a bounded escalation policy, the accession is still merged back into
one model set, and that model set participates in the existing final cohort
merge. The workflow must pass schema, lint, configuration, stub, focused
failure, and one real container-backed representative test, with exact
commands and any unverified HPC surface recorded in the handoff.
