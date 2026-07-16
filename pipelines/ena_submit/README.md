# ENA RNA-seq Alignment Evidence Submission

This pipeline submits processed RNA-seq alignments to ENA as annotation evidence.

The current production model is one ENA `ANALYSIS` per annotation/assembly partial release. Each analysis can contain many BAM/CRAM files and links back to all source runs with `RUN_REF`.

## Step 1: Build The Annotation Manifest

For an Ensembl genebuild RNA-seq directory like:

```text
.../<species>/<assembly_accession>/rnaseq
  <species>.csv
  output/<RUN>_Aligned.sortedByCoord.out.bam
```

build a manifest with:

```bash
python3 pipelines/ena_submit/bin/build_manifest_from_rnaseq_annotation.py \
  --rnaseq-dir /hps/.../annot_vert/rattus_rattus/GCA_011064425.1/rnaseq \
  --assembly-accession GCA_011064425.1 \
  --last-geneset-update 2025-12 \
  --species rattus_rattus \
  --file-format bam \
  --umbrella-study PRJEB000000 \
  --outdir /path/to/ena_manifest/GCA_011064425.1
```

Outputs:

- `manifest.tsv`: one row for the annotation-level ENA analysis.
- `files.tsv`: one row per BAM/CRAM file with run/sample/file metadata.
- `missing_files.tsv`: runs from the RNA-seq CSV without a matching alignment file.
- `summary.tsv`: counts and derived aliases.

The partial release label is derived as:

```text
<assembly_accession>-Ensembl-<last_geneset_update>
```

For example:

```text
GCA_052040795.1-Ensembl-2025-12
```

The builder currently uses the headerless RNA-seq TSV format where column 1 is sample accession, column 2 is run accession, column 9 is platform, column 11 is FASTQ URL, and column 12 is FASTQ MD5.

If `samtools` is available, the builder inspects BAM/CRAM headers for basic alignment provenance such as sort order, `@PG` programs, and STAR version. Use `--no-bam-header` to skip this.

## Step 2: Submit To ENA

Run the Nextflow workflow with the generated annotation manifest:

```bash
nextflow run pipelines/ena_submit/main.nf \
  -c pipelines/ena_submit/nextflow.config \
  --manifest /path/to/ena_manifest/GCA_011064425.1/manifest.tsv \
  --mode test \
  --webin_user "$WEBIN_USER" \
  --webin_password "$WEBIN_PASSWORD" \
  --umbrella_study PRJEB000000 \
  --remote_dir GCA_011064425.1-Ensembl-2025-12 \
  --upload_parallelism 2 \
  --outdir /path/to/ena_submit_results/GCA_011064425.1
```

Use `--mode prod` for production endpoints after Webin test validation.

### Optional BAM-to-CRAM conversion

To convert BAM inputs to CRAM before MD5 calculation and upload, enable conversion
and provide the matching reference FASTA and `.fai` index:

```bash
nextflow run pipelines/ena_submit/main.nf \
  -c pipelines/ena_submit/nextflow.config \
  --manifest /path/to/ena_manifest/GCA_011064425.1/manifest.tsv \
  --mode test \
  --convert_to_cram true \
  --reference_fasta /path/to/reference.fa \
  --webin_user "$WEBIN_USER" \
  --webin_password "$WEBIN_PASSWORD" \
  --outdir /path/to/ena_submit_results/GCA_011064425.1
```

The reference index must be at `/path/to/reference.fa.fai`. CRAM files and their
CRAI indexes are uploaded; only the CRAM files are listed in the ENA analysis XML.
Existing non-BAM inputs are passed through unchanged.

## Manifest Schema

`manifest.tsv` has one row per annotation/release:

- `files_tsv`: Path to the file-level manifest.
- `study`: Existing child study accession/alias. Leave blank to use the generated project alias.
- `project_alias`: Child project alias, usually `prj_<assembly_accession>_Ensembl_<YYYY_MM>`.
- `umbrella_study`: Umbrella study accession/alias. Recorded as analysis metadata for now.
- `analysis_alias`: Stable ENA analysis alias.
- `title`: Human-readable analysis title.
- `description`: Human-readable description.
- `assembly_accession`: INSDC assembly accession, e.g. `GCA_052040795.1`.
- `last_geneset_update`: Genome metadata value, e.g. `2025-12`.
- `partial_release_label`: Derived release label.
- `species`: Production species name.
- `taxon_id`: Optional taxon ID.
- `ref_seqs`: Optional comma-separated reference sequence accessions when assembly accession is not enough.
- `analysis_links`: Optional `Label|URL; Label2|URL2`.
- `analysis_attributes`: Optional `key=value; key2=value2`. The builder writes `attr_*` keys.
- `analysis_type`: `REFERENCE_ALIGNMENT`.
- `omit_run_refs_in_test`: `true` by default because production runs may not exist in ENA test.

`files.tsv` has one row per alignment file:

- `file_path`: Local BAM/CRAM path.
- `file_type`: `bam` or `cram`.
- `remote_name`: Filename to use in the Webin drop-box.
- `run_accession`: Source ENA/SRA/DRA run accession.
- `sample_accession`: Source sample accession from the RNA-seq CSV.
- `experiment_accession`: Optional experiment accession.
- `platform`: Sequencing platform from the RNA-seq CSV.
- `source_fastq_count`: Number of source FASTQ records for the run.
- `source_fastq_urls`: Comma-separated source FASTQ URLs.
- `source_fastq_md5s`: Comma-separated source FASTQ MD5s.
- `bam_sort_order`: Header-derived sort order when available.
- `bam_pg_programs`: Header-derived `@PG` program names when available.
- `alignment_software`: Header-derived aligner name when available.
- `alignment_software_version`: Header-derived aligner version when available.

## Outputs

- `${outdir}/ena_submission/projects/`: generated child project XML.
- `${outdir}/ena_submission/xml/`: generated analysis XML.
- `${outdir}/ena_submission/accessions.tsv`: Webin polling results.
- Nextflow work directory: upload logs and task details.

## Current Production Decisions

Locked in:

- One ENA `ANALYSIS` per annotation/assembly partial release.
- Include all source `RUN_REF` entries in production.
- Support BAM and CRAM input files.
- Derive partial release identity from assembly accession plus `last_geneset_update`.
- Use RNA-seq flat files as run/sample metadata source for now.
- Accept registry/core metadata later through CLI-compatible fields.

Pending:

- Final umbrella study accession/alias.
- Whether child project-to-umbrella linkage should be represented in ENA XML or handled externally.
- Whether production submissions should use BAMs as-is or CRAMs.
- Whether `SAMD...` sample accessions validate as `SAMPLE_REF` in Webin, or need mapping to ENA sample accessions.
- Which registry fields will replace CLI arguments for assembly/annotation metadata.
