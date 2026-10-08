# Long-read TAMA transcript models

## Runtime and module boundaries

Production runs require Nextflow 26.04.6 or newer and write timeline, report,
and trace files under `pipeline_info`. Acquisition and explicitly soft-failed
TAMA shard tasks may continue with status records; required native merge and
validation tasks terminate on failure.

The workflow is composed from these named subworkflows:

- `PREPARE_LONG_READS`: approved-manifest validation, acquisition, CCS/BAM
  conversion, and canonical FASTQ validation.
- `ALIGN_LONG_READS`: minimap2 alignment, samtools sorting/indexing, and
  alignment QC.
- `COLLAPSE_LONG_READ_MODELS` and `MERGE_LONG_READ_MODELS`: backend-native
  model generation and deterministic native/legacy merge inputs.
- `VALIDATE_COMBINED_MODELS`: canonical naming/checksum and BED12 validation.
- `RUN_DIAMOND_QC`: optional transcript, ORF, Diamond, and model-keyed report
  generation.

Model construction is selected with `--model_backend tama|stringtie2|stringtie3|tmerge|isoquant|flair|bambu|all`.
Merge semantics are selected independently with
`--backend_merge_mode native|legacy_common` (default: `native`). In native mode,
all backend processes receive the same split minimap2 BAM shards and retain
their native representation through per-accession and cohort merges:
TAMA uses BED12, StringTie2/3 use GTF, and tmerge uses native GTF. BED12 is
created only afterward for shared structural comparison. `all` produces
independent final model sets; it does not make TAMA canonical. See
`LONG_READ_BACKEND_EVALUATION.md` for the common benchmark contract and
measurements needed before selecting a backend for future Ensembl work.

`--merge_tool` is deprecated and applies only to the explicit `legacy_common`
compatibility mode. StringTie merge arguments are configured with
`--stringtie2_merge_args` and `--stringtie3_merge_args`; only the
`run_then_cohort` topology is currently supported.

Native run and cohort products emit backend-labelled model manifests, model
counts, SHA256 checksums, sorted merge file lists, and merge-file-list
checksums under `reports/native` when publishing is enabled. Set
`--publish_native_intermediates true` to publish native intermediate files as
well as the final comparison products.

The public sample contracts are `tuple val(meta), path(reads.fastq.gz)` for
canonical reads and `tuple val(meta), path(sorted.bam), path(sorted.bam.bai)`
for alignments. Helper scripts are staged as declared inputs so their exact
contents are visible in the task work directory.

This pipeline processes one complete ENA accession per alignment job, streams
minimap2 SAM output directly into `samtools sort`, collapses each accession
with TAMA, and merges models for an annotation cohort. SAM is never declared
as an output or written to disk. Sorted BAM and its index are the restartable
alignment intermediate.

The workflow has two phases. First, normalize and inspect the candidate
manifest with the classifier. If `--metadata_json` is not supplied, the
pipeline resolves ENA metadata and submitted artifacts automatically into a
persistent metadata cache, then probes the declared FASTQ headers remotely.
It preserves declared metadata separately from observed representation,
expands ENA's semicolon-separated artifacts, and writes
`run_classification.tsv`, `run_artifacts.tsv`, raw metadata snapshots, and a
cohort summary. `CONFLICT`, `UNKNOWN`, `MIXED`, and FASTQ-only raw-subread
cases are quarantined.

For PacBio, the classifier first consults NCBI Original Format metadata when
available. This is essential because ENA may expose an original file such as
`high.ccs.fq.gz` as `RUN_subreads.fastq.gz`. An Original `*.ccs.fq.gz` or
`*.ccs.fastq.gz` declaration therefore produces `PACBIO_CCS_FASTQ` and the
`ALIGN_PACBIO_CCS` direct-use action, even when normalized FASTQ headers are
archive-generated IDs. Only after submitted/original artifact evidence is
absent does the classifier use FASTQ header structure, and only then can it
reach the ungroupable raw-subread category.
Other explicit Original Format products, such as PacBio `*.flnc.fastq.gz`, are
classified as `PACBIO_PROCESSED_FASTQ` and can be aligned directly after
validation; they are not treated as raw subreads merely because ENA names its
normalized FASTQ with a `_subreads` suffix.

For new species, the candidate CSV need not be prepared manually. Supply an
NCBI taxon ID and the discovery stage will query ENA for transcriptomic ONT
and PacBio runs, enrich samples from BioSamples, retain run/sample/study and
library metadata, probe small FASTQ prefixes, group runs by BioSample, and
select a biologically diverse proposed panel. It writes a full inventory and
`proposed_long_read_manifest.tsv`; that manifest is in the versioned input
shape accepted by `normalise_manifest.py`. In taxon-only production mode, the
pipeline continues automatically after inspection and promotes only safe
`READY_FOR_REVIEW` rows into an auditable production manifest. Use `--tree` to include
subordinate taxa, `--target_sample_count` to set the soft panel size,
`--soft_download_budget` to control the default 250-GB acquisition preference,
and `--soft_raw_subread_budget` to control the default preference for one
difficult raw-subread candidate. These are soft penalties: unique biological
evidence can exceed them, with `SELECTED_WITH_WARNING` recorded explicitly.

The selector chooses one representative run per BioSample before iterative
marginal-value selection. Biological context (tissue, developmental stage,
and mixed-tissue breadth) is scored first; representation difficulty and
estimated FASTQ size are penalties, not hard exclusions. Raw PacBio subread
base counts are never treated as independent transcript coverage.

The candidate set is obtained inside the configured Ensembl genes container
with `get_transcriptomic_data.py`, then ranked by the selector above. The
default command is
`python3 -m ensembl.genes.transcriptomic_data.get_transcriptomic_data`; override
it with `--transcriptomic_data_command` if the container exposes a different
entry point.

For example:

```bash
nextflow run pipelines/long_read_tama/main.nf -profile slurm \
  --taxon_id 9913 --tree \
  --reference_fasta genome.fa \
  --fastq_cache_dir /shared/long-read-fastq-cache
```

Second, review the classification report and create an approved TSV. Production
processing accepts only rows marked `APPROVED` with `review_decision=APPROVE`;
the validator rejects unsafe paths/checksums, incompatible actions, and
unsupported classifications. The legacy whitespace parser remains a migration
adapter and no longer requires `PACBIO_SMRT`.

If human review is intentionally not required, pass `--auto_approve_safe` with
the candidate `--manifest`, `--reference_fasta`, and `--fastq_cache_dir`. The
pipeline then automatically promotes every `READY_FOR_REVIEW` row with a
compatible classification, quarantines the remaining rows, writes
`automatic_selection_audit.tsv`, and continues with production processing. It
fails if no safely runnable rows are found. This mode is opt-in; without it,
`--manifest` remains inventory-only.

`fastq-dl=2.0.1` downloads one accession at a time into a persistent directory
provided with `--fastq_cache_dir`. A lock prevents concurrent writers for the
same filename. The resulting single FASTQ is gzip-checked, compared with the
manifest checksum, and checked against the manifest filename before an atomic
rename into the cache. A valid cached file is reused on later runs; an invalid
cached file fails rather than being silently overwritten.

By default, sorted BAMs are partitioned by reference contig before TAMA
Collapse. TAMA runs per contig, and those beds are merged per accession before
the cohort merge. Use `--shard_mode none` only for deliberately small inputs.
Arbitrary read sharding and genomic windows are intentionally unsupported.

If an individual TAMA shard has been deliberately abandoned, it can be omitted
without dropping the rest of its accession by passing exact comma-separated
`run:shard` keys, for example
`--skip_tama_shards SRR32588732:10,SRR32588732:MT`. The remaining successful
shards are still included in the accession and cohort merges.

Known per-shard TAMA failures are recorded in `tama_status.tsv` and do not
produce a BED for that shard; successful shards continue to the accession and
cohort merges. TAMA exit statuses are captured by the shard adapter and are
non-terminal. The separate `--merge_tool` option controls only the legacy
common post-collapse accession/cohort merge stage; it is not the
model-construction backend selector. For a direct native backend comparison,
use `--model_backend all --backend_merge_mode native`.

The TAMA adapter has a targeted recovery for the known TAMA 1.0.3 empty-locus
crash (`IndexError: list index out of range`). It retries only that signature
after filtering the shard to primary mapped alignments (`samtools -F 2308`),
records the retry in `tama_collapse.stderr`, and still records a non-zero
status if recovery fails. This is a compatibility workaround, not a biological
equivalence guarantee: supplementary/secondary evidence is removed on retry
and recovered model counts should be reviewed.

`tmerge` consumes read-level exon GTF records. The `model_backend=tmerge` path
therefore converts each split minimap2 BAM directly to read-level GTF before
calling tmerge; it does not convert already-collapsed TAMA BED12 models. It is
not a drop-in replacement for TAMA's collapse algorithm.

TAMA allocations use configurable workload tiers. The defaults are
`128.GB`, `256.GB`, and `512.GB`. The TAMA adapter records exit statuses such
as 137, 140, and 143 in `tama_status.tsv`; other
tool and validation failures are also non-terminal at the process-policy
boundary and remain visible in the trace/report outputs. `--shard_contig_reads` is an operational mapped-read
heuristic, not a biological or guaranteed memory boundary.

Example:

```bash
nextflow run pipelines/long_read_tama/main.nf -profile slurm \
  --approved_manifest approved_run_manifest.tsv \
  --fastq_cache_dir /shared/long-read-fastq-cache \
  --reference_fasta genome_dumps/bos_indicus_toplevel.fa \
  --outdir /path/to/results
```

The inventory tool can be run without acquiring complete reads when metadata
and selected artifacts are represented by local/static fixtures:

```bash
python3 pipelines/long_read_tama/bin/read_input_classification.py inspect \
  normalised_manifest.tsv metadata.json classification_report
```

For live preflight, ENA responses are cached losslessly per accession and
reused on later inventory runs. SRA evidence is recorded as `unavailable` when
the toolkit is absent. `--metadata_json` can still be used to replay a frozen
inventory without network access.

```bash
python3 pipelines/long_read_tama/bin/resolve_metadata.py \
  accessions.txt /shared/long-read-metadata-cache metadata.json
```

Run a structural test with:

```bash
nextflow run pipelines/long_read_tama/main.nf -stub-run -profile stub \
  --approved_manifest pipelines/long_read_tama/test/approved_run_manifest.tsv \
  --fastq_cache_dir /tmp/long-read-tama-cache \
  --reference_fasta pipelines/long_read_tama/test/reference.fa
```

An approved structural run can use `pipelines/long_read_tama/test/approved_run_manifest.tsv`
with `--approved_manifest`; ONT rows use `splice`, while PacBio CCS rows use
`splice:hq`. PacBio BAM routes require a site-validated CCS tool, sidecars,
arguments, and canonical FASTQ contract. Provide `--ccs_container`,
`--ccs_command`, `--ccs_args`, and `--ccs_expected_version`; the BAM process
verifies the version before converting and records it in `molecule_audit.tsv`
and `versions.yml`.

The pipeline pins GS-TAMA 1.0.3 by default. Set `--tama_container` only when
using a site-local mirror or validated replacement image. The `splice:hq` /
secondary-alignment choice must still be benchmarked on
`SRR29278220_subreads.fastq` before rollout.

## IsoQuant backend

IsoQuant is an independent backend over the complete sorted/indexed BAMs. It
does not consume contig shards and its native GTF, read evidence, counts, and
audit manifest are kept separate from TAMA/StringTie/tmerge merging. The
default IsoQuant mode is annotation-free cohort discovery:

```bash
nextflow run pipelines/long_read_tama/main.nf -profile slurm \
  --approved_manifest approved_run_manifest.tsv \
  --fastq_cache_dir /shared/long-read-fastq-cache \
  --reference_fasta genome.fa \
  --model_backend isoquant \
  --isoquant_data_type pacbio_ccs \
  --isoquant_scope cohort
```

Use `--isoquant_scope accession` or `both` for per-accession evidence. The
reference-guided mode requires `--isoquant_genedb`; annotation-free mode
rejects that parameter. The pinned BioContainers image is IsoQuant 4.0.0 and
the process uses the upstream-supported `isoquant --reference --bam
--data_type --analysis --large_output` command shape.
Under the local Apptainer runtime, the workflow redirects IsoQuant's `HOME`
to a writable task-local directory when the image home is read-only. The
optional `--isoquant_numba_disable_jit true` flag is retained only as a
diagnostic fallback and is not the normal performance setting.
In annotation-free mode, IsoQuant's native `transcript_model_reads.tsv.gz` is
used as the stable read-evidence fallback when a separate `read_info.tsv.gz`
is not emitted; the audit manifest records that source explicitly.

## FLAIR backend

FLAIR is an independent annotation-free transcript-discovery backend. It first
derives a junction BED from each complete accession BAM with the repository's
CIGAR-aware BAM helper, then runs `flair transcriptome` and combines the native
accession products with `flair_combine`. The helper retains splice `N`
operations while accepting the clipping and indel operations produced by
minimap2. No annotation is supplied. Stable products include the native
BED/GTF/FASTA and read-to-isoform map; the cohort BED12 is derived only after
FLAIR's native combine stage. Set
`--model_backend flair` or include it in `--model_backend all`.

FLAIR and IsoQuant intentionally use complete sorted/indexed BAMs rather than
contig shards. They must be compared with the same BAMs and explicit library
technology as the other backends, but should not be treated as having TAMA or
StringTie merge semantics.

## Bambu backend

Bambu is included as an additional annotation-free complete-BAM discovery
backend. It runs with `annotations=NULL` and `NDR=1`; quantification and
read-tracking are enabled by default so the candidate models retain native
read evidence. Set `--model_backend bambu` or include it in
`--model_backend all`; use `--bambu_scope both` for accession and cohort
outputs. The pinned Bioconda image omits `BiocManager`, so the default
`--bambu_install_biocmanager true` bootstraps that small dependency into the
task-local R library. On an executor without outbound package access, provide
a site-built `--bambu_container` that already contains BiocManager and set
the bootstrap flag to false. Alternatively, install the dependency once in a
shared library and pass its path with `--bambu_r_lib`, then set
`--bambu_install_biocmanager false`.
Read-to-transcript assignments are retained by default with
`--bambu_track_reads true`; disable this only for a deliberately
structural-only, lower-memory comparison. Use `--bambu_quantify false` for
discovery-only output, in which case read assignments are not expected.

## Weekend comparison run

For an ONT cohort, the primary exploration run is:

```bash
nextflow run pipelines/long_read_tama/main.nf -profile slurm \
  -with-trace -with-report -with-timeline \
  --approved_manifest approved_run_manifest.tsv \
  --fastq_cache_dir /shared/long-read-fastq-cache \
  --reference_fasta genome.fa \
  --model_backend all \
  --backend_merge_mode native \
  --backend_failure_policy continue \
  --isoquant_data_type nanopore \
  --isoquant_scope both \
  --bambu_scope both \
  --outdir /shared/long-read-results/all-native
```

For PacBio CCS, change only `--isoquant_data_type` to `pacbio_ccs` and ensure
the approved manifest classifies the reads accordingly. The primary cohort
and accession candidate sets are published under backend-specific directories;
native tool products and evidence are retained alongside them. The
backend-neutral structural summaries are published as
`reports/comparison/candidate_model_comparison_cohort.tsv/.json`. Accession
summaries are run independently per accession and published as
`candidate_model_comparison_accession_<accession>.tsv/.json`, with matching
`candidate_model_manifest_accession_<accession>.tsv` files. `backend_failure_policy=continue`
allows remaining methods to finish if an optional backend fails; the failure
remains visible in the trace and must be included in the Monday review.
Summarise computational cost after the run with:

```bash
python3 pipelines/long_read_tama/bin/summarise_backend_resources.py \
  /shared/long-read-results/all-native/pipeline_info/execution_trace.txt \
  /shared/long-read-results/all-native/reports/comparison/backend_resources.tsv \
  /shared/long-read-results/all-native/reports/comparison/backend_resources.json
```

Do not compare an ONT run with the default IsoQuant `pacbio_ccs` setting. The
technology is intentionally explicit so that a filename or accession label
cannot silently change the model-building behaviour.

## Combined-model Diamond validation

The optional combined-model validation path starts after each backend's
cohort-wide native merge:

```text
combined_models.bed -> combined_transcripts.fa -> combined_transcripts.faa -> Diamond
```

Enable it with an existing Diamond database:

```bash
nextflow run pipelines/long_read_tama/main.nf -profile slurm \
  --approved_manifest approved_run_manifest.tsv \
  --fastq_cache_dir /shared/long-read-fastq-cache \
  --reference_fasta genome.fa \
  --run_diamond_validation \
  --diamond_reference_db reference.dmnd
```

Alternatively, pass `--diamond_reference_proteins` and the pipeline will build
the database with the pinned Diamond module. The initial peptide strategy is
the longest ATG-initiated ORF per transcript. It emits
`combined_orf_manifest.tsv`, including `NO_ATG` and `PARTIAL_ORF` rows, so a
future ORF predictor can replace `PREDICT_LONGEST_ATG_ORFS` while preserving
the same model-keyed Diamond report interface. Diamond reports retain every
combined model and assign `HIT`, `NO_PROTEIN_HIT`, or `NO_PEPTIDE` status.

The report also emulates the legacy `HiveClassifyTranscriptSupport` coverage
and identity classes. The default `--diamond_classification_type long_read`
uses classes `1`–`7` with the historical coverage/identity thresholds and
emits the corresponding `_1`–`_7` suffix in
`legacy_biotype_suffix`. The original Perl runnable updates an Ensembl
transcript biotype in SQL; this pipeline records the suffix as a report field
because the combined TAMA models are not yet in an Ensembl core database.
