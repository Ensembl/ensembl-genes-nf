# Funannotate CRAM-backed prediction

This pipeline runs one `funannotate predict` analysis per tab-separated
manifest row. The genome and all annotation/evidence files are local files on
the HPS filesystem. Only the RNA-seq CRAMs are downloaded from FTP. The
pipeline stages the inputs, prepares derived evidence, converts CRAMs to BAM,
and runs Funannotate in a container.

The default executor is Slurm. A `local` profile is available for development
and stub checks.

## Workflow

The orchestration is implemented in `subworkflows/run_funannotate_pipeline.nf`:

1. Parse and validate the manifest.
2. Stage the local genome and optional local evidence files.
3. Download the CRAMs from FTP and verify requested SHA-256 checksums.
4. Create a genome FASTA index with `samtools faidx`.
5. Check every CRAM with `samtools quickcheck`.
6. Convert each CRAM to BAM using the staged genome as the reference.
7. Merge the converted BAMs, sort the merged BAM, and create its index. The
   main BAM passed to Funannotate is `staged/all_runs.merged.sorted.bam`.
8. Convert all protein-alignment GTFs into one Funannotate GFF3 file.
9. Convert the Scallop GTF to transcript FASTA with `gffread` when present.
10. Run `funannotate predict` using the merged BAM and available evidence.

The original staged inputs are retained. Derived files are added to the staged
directory for the Funannotate task.

## Manifest

The manifest is tab-separated and must contain these columns:

* `cohort_id`
* `assembly_accession`
* `species`
* `genome_path`
* either `cram_ftp_urls` or `cram_ftp_dir`

Optional columns are:

* `cram_sha256s`
* `cram_ftp_dir`
* `transcript_evidence_paths`
* `stringtie_gtf_path`
* `scallop_gtf_path`
* `protein_alignment_paths`
* `trnascan_path`

`cram_ftp_urls` and all `*_paths` lists are semicolon-separated. CRAM URLs may
use either `ftp://` or `https://`. Alternatively, `cram_ftp_dir` may contain a
directory URL; the pipeline then selects only files matching
`*.Aligned.SortedByCoordinate.out.cram`. CRAM checksums
must be in the same order as `cram_ftp_urls` when supplied. Local paths may be
absolute HPS paths or paths relative to the directory from which Nextflow is
launched. A gzip-compressed genome FASTA (`.fa.gz`) is decompressed into the
staging area as `genome.fa` because `samtools faidx` requires an uncompressed
or BGZF-compressed FASTA.

`protein_alignment_paths` may contain multiple local GTF paths, for example:

```text
/hps/.../genblast1.gtf;/hps/.../genblast2.gtf
```

The source basenames are retained in generated alignment IDs. Each input GTF
must contain `transcript` and `exon` records with a `transcript_id` attribute.

## Evidence handling

### StringTie and Scallop

The StringTie GTF is staged as `stringtie.gtf` and passed directly using
`--stringtie`.

The Scallop GTF is staged as `scallop.gtf`. When present, `gffread` uses the
staged genome to create `staged/scallop_transcripts.fa`. This FASTA is passed as
`--transcript_evidence`.

### Protein alignments

Each protein-alignment GTF is staged separately. The script
`bin/convert_protein_alignments.py` reproduces the manual GenBlast conversion
and writes `staged/genblast_protein_alignments.gff3`. That combined GFF3 is
passed using `--protein_alignments`.

### RNA-seq CRAMs

Each FTP CRAM is converted to a BAM with `samtools view -T genome.fa`. The BAMs
are then merged, sorted, and indexed. Funannotate receives one BAM argument:
`staged/all_runs.merged.sorted.bam`. Set `--cpus` to the requested thread
count for CRAM conversion, merging, sorting, and Funannotate.

## Containers

Containers are declared in the individual modules, following the repository's
module convention. Nextflow downloads them through the configured Singularity
runtime; no local absolute image paths are required.

| Process | Container |
|---|---|
| FTP staging | `quay.io/biocontainers/curl:7.80.0` |
| CRAM validation, conversion, and BAM merging | `quay.io/biocontainers/samtools:1.20` |
| Scallop GTF conversion | `docker://getwilds/gffread:0.12.7` |
| Protein-alignment conversion | Repository-standard `python` process container |
| Funannotate prediction | `docker://nextgenusfs/funannotate` |

The local profile disables containers because stub mode does not execute the
real commands.

## Running

### Slurm

Slurm is the default executor. The explicit profile is also available:

```bash
nextflow run pipelines/funannotate/main.nf \
  -profile slurm \
  --input_manifest /path/to/funannotate.tsv \
  --stage_dir /hps/nobackup/flicek/ensembl/genebuild/user/funannotate-stage \
  --outdir /path/to/results \
  --cpus 8 \
  -resume
```

The pipeline supplies only the Funannotate arguments listed in the workflow
description: input genome, output directory, species, merged RNA-seq BAM,
StringTie evidence, transcript evidence, protein alignments, tRNAscan, and
CPUs. The Funannotate database is provided by the container.

### Local stub check

The local profile uses the local executor, disables containers, and uses stub
implementations for all processes:

```bash
nextflow run pipelines/funannotate/main.nf \
  -stub-run \
  -profile local \
  --input_manifest pipelines/funannotate/test/data/funannotate_manifest.tsv \
  --stage_dir /tmp/funannotate-stage \
  --outdir /tmp/funannotate-stub
```

## Outputs

Results are published under:

```text
<outdir>/<assembly_accession>/<safe_cohort_id>/
```

The main output is the Funannotate prediction directory, including the
resulting GFF3 annotation and Funannotate logs. Each process produces a
`versions.yml` file for tool-version tracking. FTP inputs and checksums are
retained under the configured `stage_dir`. The local source files are not
modified.
