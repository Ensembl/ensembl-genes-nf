# Long-read TAMA transcript models

## Runtime and module boundaries

Production runs require Nextflow 26.04.6 or newer and write timeline, report,
and trace files under `pipeline_info`. Required stages terminate the run on
failure; retries are limited to network/cache acquisition.

The workflow is composed from these named subworkflows:

- `PREPARE_LONG_READS`: approved-manifest validation, acquisition, CCS/BAM
  conversion, and canonical FASTQ validation.
- `ALIGN_LONG_READS`: minimap2 alignment, samtools sorting/indexing, and
  alignment QC.
- `COLLAPSE_LONG_READ_MODELS` and `MERGE_LONG_READ_MODELS`: TAMA model
  generation and deterministic merge inputs.
- `VALIDATE_COMBINED_MODELS`: canonical naming/checksum and BED12 validation.
- `RUN_DIAMOND_QC`: optional transcript, ORF, Diamond, and model-keyed report
  generation.

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

Known per-shard TAMA failures are recorded in `tama_status.tsv` and do not
produce a BED for that shard; successful shards continue to the accession and
cohort merges. Missing executables and scheduler/resource kills remain fatal or
retryable. The default merge backend is TAMA. An experimental `tmerge` backend
can be selected with `--merge_tool tmerge`; it runs as a separate `TMERGE`
process using the pinned public
`community.wave.seqera.io/library/pip_tmerge:6cf60ff0bf166552` image. `tmerge`
consumes a coordinate-sorted exon GTF, so the pipeline converts collapsed BED12
models to GTF and converts its output back to BED12. In contig-sharded mode the
selected backend is also used for the per-accession regrouping step before the
cohort merge. It is not a drop-in replacement for TAMA's collapse algorithm.

TAMA allocations use configurable workload tiers. The defaults are
`128.GB`, `256.GB`, and `512.GB`, with five bounded retries. Exit statuses
137, 140, and 143 advance through the tiers and are ignored after the retry
budget is exhausted; validation, malformed-input, and other tool failures
still terminate. `--shard_contig_reads` is an operational mapped-read
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

## Combined-model Diamond validation

The optional combined-model validation path starts only after the cohort-wide
TAMA merge:

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
