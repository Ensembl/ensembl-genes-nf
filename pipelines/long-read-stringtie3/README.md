# Long-read alignment and de novo StringTie3

This focused pipeline aligns long-read RNA-seq FASTQs to selected genome
assemblies, retains each sorted/indexed BAM, performs de novo StringTie3
assembly per run, merges transcript models per species/assembly, then estimates
abundance for each run against that StringTie-derived merged catalog. It does
not consume external gene annotations.

StringTie performs the transcript assembly, merge, and quantification. Workflow
code handles input validation, reference indexing, alignment, and scheduling.

## Run manifest

Provide a tab-delimited file with one row per run:

```tsv
species_id	assembly_id	assembly_accession	run_id	sample_id	technology	fastq	minimap2_preset	nascent_mode
abramis_brama	GCA_963993115.1	GCA_963993115.1	ERR13757885	SAMEA11296570	pacbio	/data/fastq/ERR13757885.fastq.gz	splice:hq	off
```

Required columns are `species_id`, `assembly_id`, `assembly_accession`,
`run_id`, `technology`, and `fastq`. `sample_id`, `minimap2_preset`, and
`nascent_mode` are optional. FASTQ may be a local path (absolute or relative to
the manifest) or a remote HTTP(S)/FTP URI. The bare `ftp.sra.ebi.ac.uk/...`
form in project batch manifests is interpreted as HTTPS. The workflow fetches
one genome per species/assembly using `bin/fetch_genome.py`: it prefers the
latest available Ensembl Organisms release, falls back to NCBI Datasets, and
normalizes FASTA headers using the selected source's rules before indexing.
The NCBI assembly report is retained beside the fetched FASTA as provenance.

A FASTQ column identifies one input file per run. For multiple lanes from one
biological sample, supply one row per run and use `sample_id` as the shared
metadata key; outputs remain per run. `run_id` must be unique across the
manifest. IDs may contain letters, numbers, dots, underscores, and hyphens.
Technology is `ont` or `pacbio`; default minimap2 presets are `splice` and
`splice:hq`, respectively. The supported presets share splice index settings.

All rows for one `species_id` + `assembly_id` must use the same
`assembly_accession` and nascent mode. An assembly gets one genome download and
one minimap2 index, shared by its runs. The downloader normalizes reference
sequence names before alignment and transcript assembly.

`nascent_mode` defaults to `off`. For total-RNA / rRNA-depleted libraries,
StringTie3 recommends nascent-aware assembly: set it to `mature` to use `-N`
(emit mature models only), or `retain` to use `--nasc` (also retain nascent
intermediates in the GTF). Keep the value consistent within an assembly so runs
contribute comparable transcript sets to the merge.

## Outputs

```text
pipeline_info/
  input_manifest.tsv
  input_manifest.tsv.sha256
per_assembly/<species>/<assembly>/reference/
  GCA_...vN.fa
  *_assembly_report.txt (when supplied by NCBI)
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
per_run/<species>/<assembly>/<run>/stringtie3/quantification/
  quantified.gtf
  gene_abundances.tsv
  ballgown/*.ctab
```

The first StringTie pass uses `-L` to assemble long reads de novo and `-A` to
write a gene abundance table. `--merge` combines each assembly's per-run GTFs
without adding external annotation. The second pass uses `-L -e -G` with that
StringTie-generated merged GTF, so every run is quantified against the same
transcript catalog; `-A` and `-b` produce gene abundance and Ballgown coverage
tables. The shared-catalog abundance table is the one to use for cross-run
comparisons. StringTie creates these GTF and evidence tables directly.

The retained BAM and CSI are copied to the output directory while a work copy
remains available for downstream tasks and resume. Genome downloads are
published once per assembly; each assembly's fetch and index can proceed while
other assemblies are still downloading. Budget scratch for BAMs across the run,
their published copies, staged references, and temporary sort files. Set
`-work-dir` to scratch. After confirming published outputs,
`nextflow clean` can reclaim the work copy; cleaning it removes resume data.

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
`-profile docker_runtime` for Docker, or `-profile apptainer` where Apptainer
is installed. Defaults are two concurrent genome downloads, index builds,
alignments, StringTie assemblies, and quantification tasks, with one merge at a
time. Set `max_genome_download_jobs`, `max_index_jobs`, `max_alignment_jobs`,
`max_stringtie_assembly_jobs`,
`max_stringtie_merge_jobs`, and `max_stringtie_quantification_jobs` to match
the site allocation and scratch capacity. The per-stage caps bound large tasks
in flight; `executor_queue_size` provides an additional total cap on queued or
running tasks for grid executors.

For a wiring-only check with the bundled two-species fixture:

```bash
nextflow run pipelines/long-read-stringtie3/main.nf -stub-run -profile test \
  --manifest pipelines/long-read-stringtie3/test/data/multi_species.tsv \
  --outdir /tmp/long-read-stringtie3-test
```

Stub execution checks channel wiring and output layout; use a real, small
container-backed run before production deployment.
