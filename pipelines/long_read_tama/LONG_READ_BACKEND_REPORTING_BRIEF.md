# Long-read backend results reporting brief

## Purpose

Produce a reusable report comparing the implemented long-read model backends and their merge/QC behaviour. The report must work for the current *Danio rerio* runs and for later, more complete reruns without treating a partially populated output directory as a completed experiment.

The report should answer two separate questions:

1. **What did each backend produce?**
2. **How trustworthy and comparable are those products for this particular run?**

## Current inventory snapshot

The HPC output trees currently provide the following useful material:

- Nine accession-level sorted/indexed alignment BAMs, alignment statistics, FASTQ provenance, and validated manifests.
- Extensive per-shard native model output for StringTie2, StringTie3, tmerge, and TAMA, including GTF/BED products and shard status files.
- IsoQuant native products, discovered transcript/gene counts and TPMs, transcript models, comparison BED12, and status/manifests.
- Bambu native products including models, read assignments, counts, extended annotations, and the Bambu RDS object.
- FLAIR junction, transcript, and BED12 products in the published tree, but older FLAIR products must be marked provisional until regenerated after the contig-name and GNU-sort fixes.
- Accession/cohort merge products, merged BED12 files, merge filelists, merged gene/transcript reports, and model validation reports in the native/skip-TAMA trees.
- A separate Diamond/QC tree containing backend-specific Diamond tables for StringTie2, StringTie3, TAMA, and tmerge.
- Nextflow execution reports, timelines, traces, software versions, validation JSON/TSV files, shard plans, and backend status tables.

The inventory is sufficient to report alignment availability, backend coverage, shard completion, model counts, merge behaviour, and resource usage. It is not yet sufficient to claim that every backend is directly comparable: active TAMA shards, ignored/failed tasks, stale outputs from earlier revisions, and the pre-fix FLAIR products need to be separated first.

## Provenance rules

Every report must record, for each result tree:

- results directory and HPC host/path;
- Nextflow version, run name, launch/completion time, profile, and repository revision;
- manifest, reference FASTA, minimap2 index, and their checksums where available;
- complete parameter set, including backend scopes, merge mode, failure policy, and skip lists;
- container/image identifiers and software versions;
- whether the result tree is complete, still running, resumed, cached, or assembled from more than one run.

Do not infer completion from the presence of a BED/GTF file alone. Published directories can retain products from previous attempts. A product is reportable only when its corresponding process/shard is `SUCCESS` or `CACHED` in the execution trace/status records for the same run and revision.

Use these status categories consistently:

`SUCCESS`, `CACHED`, `RUNNING`, `FAILED`, `IGNORED`, `SKIPPED`, `EMPTY_OUTPUT`, `STALE_OR_UNVERIFIED`.

`EMPTY_OUTPUT` is different from `FAILED`: an executed tool can legitimately produce no models, but that must be demonstrated by the tool’s status/QC output rather than assumed from a zero-byte or absent file.

## Required report sections

### 1. Run and input summary

Summarise accessions, read source/platform, read counts, mapped/primary/supplementary reads, mapping rate, reference, and shard plan. Include the number of contigs/shards and the exact skipped shard keys. Show a compact run matrix so the reader can see which accessions and scopes are complete.

### 2. Execution/completeness matrix

Create a matrix with rows for backend, scope, accession/cohort, and shard, and columns for status, output path, process/task ID, attempt, runtime, CPU, memory, and error/ignore reason. This is the authoritative answer to “what is available”.

At minimum include:

- alignment and alignment QC;
- TAMA collapse and merge;
- StringTie2 collapse and merge;
- StringTie3 collapse and merge;
- tmerge collapse and merge;
- IsoQuant accession/cohort products;
- Bambu accession/cohort products;
- FLAIR junctions, transcriptome, GTF-to-BED12, and combine;
- canonicalisation, validation, software-version collection, and optional Diamond QC.

### 3. Backend output summary

For every backend and scope, report:

- model/transcript/gene counts;
- exon and splice-junction counts where available;
- read-assignment/count products where available;
- per-accession and cohort products separately;
- native output format and normalised BED12 output;
- number of successful, failed, ignored, skipped, and empty shards;
- tool version, container, key parameters, wall time, CPU hours, maximum RSS, and cache status.

Do not compare raw model counts without stating whether they are per-shard, per-accession, cohort-merged, or post-canonicalisation counts.

### 4. Merge topology and model comparison

Explain the actual path taken through the workflow:

```text
reads -> alignment -> contig shards -> backend collapse
      -> accession merge -> cohort merge -> canonicalise/validate
      -> backend comparison and optional Diamond QC
```

Make clear that StringTie2, StringTie3, tmerge, TAMA, IsoQuant, Bambu, and FLAIR do not necessarily emit equivalent intermediate objects. Compare at the common normalised BED12/canonical model layer where possible, while retaining native outputs for tool-specific interpretation.

For each pair or multi-way comparison, report shared, backend-specific, and unresolved models; splice-junction agreement; transcript-length/exon-number distributions; and the effect of accession versus cohort merging. Include the denominator and filtering rules for every overlap percentage.

### 5. QC and failure interpretation

Summarise validation results, empty outputs, ignored errors, failed tasks, missing shards, and known technical limitations. Include the relevant `.command.err`/trace reference for each issue.

Known run-specific cautions:

- FLAIR results produced before the corrected Wave container and BAM-contig-preserving junction extraction must be labelled provisional and excluded from final biological comparisons until rerun.
- Bambu cohort output from the old parallel execution must be labelled failed/ignored and rerun with the corrected single-core reducer workaround; successful accession-level Bambu products remain reportable.
- TAMA counts and comparisons must remain provisional while any TAMA collapse shards are still running or explicitly skipped.
- A directory named `skip_tama` must not be interpreted from its name alone; verify the execution trace and process statuses because cached/stale TAMA products may still be present.

### 6. Resource and performance summary

Use the Nextflow trace and timeline to show runtime, CPU hours, maximum RSS, queue/wait time if available, retries, cache rate, and task counts by process/backend. Rank the slowest shards and state whether slowness is caused by read volume, contig size, tool behaviour, queue time, or downstream serialisation/merge.

This section should distinguish active work from completed work and should never count a running shard as a zero-result backend.

### 7. Biological/tool-use interpretation

End with a restrained “what each tool is good for” section based on observed outputs and failure modes, not just model totals. For example:

- alignment/QC: input and mapping quality;
- StringTie2/StringTie3: reference-guided or assembly-style transcript reconstruction and fast shard-level comparison;
- tmerge: merging transcript models across accessions/shards;
- TAMA: long-read-oriented collapse/merge behaviour, with runtime and shard sensitivity explicitly reported;
- IsoQuant: native abundance and discovered-model products;
- Bambu: model and read-assignment/count products, noting cohort parallelism limitations in the affected version;
- FLAIR: long-read transcript reconstruction, only after the corrected container and contig-compatible junction input are confirmed;
- Diamond: optional protein/homology support for model QC, not a replacement for transcript-structure comparison.

## Recommended machine-readable deliverables

The report-generating session should emit these stable tables, using TSV/CSV plus a human-readable Markdown or HTML summary:

- `run_manifest.tsv` — one row per run/result tree and provenance fields;
- `input_alignment_summary.tsv` — one row per accession;
- `backend_status.tsv` — one row per backend/scope/accession/cohort;
- `shard_status.tsv` — one row per backend/accession/shard;
- `backend_model_counts.tsv` — counts with level (`shard`, `accession`, `cohort`, `canonical`);
- `comparison_summary.tsv` — pairwise/multi-way overlap and junction comparison metrics;
- `merge_summary.tsv` — merge inputs, outputs, counts, and reports;
- `failure_summary.tsv` — ignored/failed/empty/stale products and evidence paths;
- `resource_summary.tsv` — runtime/resource/cache statistics;
- `software_versions.tsv` — process/container/tool versions.

Each table should contain `run_id`, `repo_revision`, `result_root`, and `status` where applicable. Preserve unknown or unavailable values as `NA`; do not silently coerce them to zero.

## Reuse and acceptance criteria

The future reporting session should accept one or more result roots, for example:

```text
--results-root results/danio_rerio_all_native
--results-root results/danio_rerio_all_native_skip_tama
--results-root results/danio_rerio_with_diamond
```

It should compare only common accession/scope/backend units unless the report explicitly labels a comparison as incomplete. It should be safe to rerun after the remaining jobs finish, and it should produce the same table schema for a small test run, the current partial run, and a complete future run.

The report is accepted when a reader can determine, without inspecting the work directory manually:

1. which backend/scope/shards completed;
2. which outputs are valid, provisional, stale, empty, failed, ignored, or skipped;
3. how many models and junctions each backend produced at each merge level;
4. how the backends agree and differ on a common normalised representation;
5. which results are affected by known FLAIR/Bambu/TAMA issues;
6. which backend is slow or resource-intensive and why; and
7. exactly which commands, revisions, containers, and input files produced the numbers.

