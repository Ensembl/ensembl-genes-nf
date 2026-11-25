# Long-read alignment and StringTie3

This focused pipeline aligns long-read RNA-seq FASTQs to multiple species or
assemblies in one Nextflow execution, retains each sorted/indexed BAM, runs
StringTie3 per run, merges transcript models per assembly, then quantifies each
run against that assembly's merged transcript set.

It uses StringTie for transcript assembly, merging, gene abundance reporting,
fully covered reference transcript reporting, and abundance/coverage tables.
There are no custom transcript assembly, merge, or quantification scripts.

## Run manifest

Provide a tab-delimited file with one row per run:

```tsv
species_id	assembly_id	run_id	sample_id	technology	fastq	reference_fasta	annotation_gtf	minimap2_preset
Homo_sapiens	GRCh38	run001	sample001	pacbio	reads/run001.fastq.gz	refs/GRCh38.fa	refs/GRCh38.gtf	splice:hq
Mus_musculus	GRCm39	run002	sample002	ont	reads/run002.fastq.gz	refs/GRCm39.fa	refs/GRCm39.gtf	splice
```

Required columns are `species_id`, `assembly_id`, `run_id`, `technology`,
`fastq`, `reference_fasta`, and `annotation_gtf`. `sample_id` and
`minimap2_preset` are optional. Paths may be absolute or relative to the
manifest. A FASTQ column identifies one input file per run. For multiple lanes
from one biological sample, supply one row per run and use `sample_id` as the
shared metadata key; the workflow still produces per-run outputs. IDs may
contain letters, numbers, dots, underscores, and hyphens. Technology is `ont`
or `pacbio`; the corresponding default presets are `splice` and `splice:hq`.

All rows for one `species_id` + `assembly_id` must use the same reference FASTA
and annotation. An assembly gets one minimap2 index, shared by its runs. The
index is built with splice-mode settings; the supported `splice` and
`splice:hq` mapping presets share those index settings. Reference annotation
is required so StringTie can produce guided assemblies, fully covered
reference transcript GTFs, and quantification tables.

## Outputs

```text
pipeline_info/
  input_manifest.tsv
  input_manifest.tsv.sha256
per_assembly/<species>/<assembly>/reference/
  reference.mmi
  reference.fa.fai
per_assembly/<species>/<assembly>/stringtie3/merge/
  merged.gtf
per_run/<species>/<assembly>/<run>/alignment/
  aligned.sorted.bam
  aligned.sorted.bam.csi
  flagstat.txt
  idxstats.tsv
per_run/<species>/<assembly>/<run>/stringtie3/assembly/
  <run_id>.assembled.gtf
  gene_abundances.tsv
  fully_covered_reference.gtf
per_run/<species>/<assembly>/<run>/stringtie3/quantification/
  quantified.gtf
  gene_abundances.tsv
  ballgown/*.ctab
```

The first StringTie pass uses `-L -G -A -C` to discover models and report gene
abundance and fully covered reference transcripts. `--merge -G` produces the
per-assembly unified transcript set. The second pass uses `-L -e -G -A -b` to
estimate abundance for every run on that shared set and emits Ballgown coverage
tables directly from StringTie. The two abundance tables answer different
questions: discovery-pass abundance is tied to each run's own assembled set;
quantification-pass abundance uses the common merged set and is the one to use
for cross-run comparisons.

The retained BAM and CSI are copied to the output directory while a work copy
remains available for downstream tasks and resume. Budget scratch for the BAMs
across the whole run, those published copies, and temporary sort files. Set
`-work-dir` to scratch. After confirming published outputs, `nextflow clean`
can reclaim the work copy; cleaning it removes resume data.

## Run

```bash
nextflow run pipelines/long-read-stringtie3/main.nf \
  --manifest /data/long_reads.tsv \
  --outdir /data/long-read-results \
  -work-dir /scratch/nextflow/long-read-stringtie3 \
  -profile singularity \
  -with-report /data/long-read-results/pipeline_info/report.html \
  -with-trace /data/long-read-results/pipeline_info/trace.tsv
```

The repository's shared configuration defaults to Singularity. Use
`-profile docker_runtime` for Docker, or `-profile apptainer` where Apptainer is
installed. Defaults are two concurrent index builds, alignments, StringTie
assemblies, and quantification tasks, with one merge at a time. Set
`max_index_jobs`, `max_alignment_jobs`,
`max_stringtie_assembly_jobs`, `max_stringtie_merge_jobs`, and
`max_stringtie_quantification_jobs` to match the site allocation and scratch
capacity. These per-stage caps bound the number of large BAM and StringTie
tasks in flight while allowing work for independent assemblies to overlap.
`executor_queue_size` provides an additional total cap on queued/running tasks
for grid executors; its default is eight. The local executor instead observes
available CPU and memory alongside the per-stage caps.

For a wiring-only check with the bundled two-species fixture:

```bash
nextflow run pipelines/long-read-stringtie3/main.nf -stub-run -profile test \
  --manifest pipelines/long-read-stringtie3/test/data/multi_species.tsv \
  --outdir /tmp/long-read-stringtie3-test
```

Stub execution checks channel wiring and output layout; use a real, small
container-backed run before production deployment.
