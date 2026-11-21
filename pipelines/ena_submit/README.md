# ENA RNA-seq Alignment Evidence Submission

This pipeline submits processed RNA-seq alignments to ENA as annotation evidence.

## What this pipeline does

An Ensembl genebuild produces many alignment files, but ENA submission is
organised around ENA projects, samples, runs, and analyses. This pipeline is
the translation layer between those two models. It takes one annotation-level
manifest plus a file-level manifest, validates and prepares the alignments,
then creates and submits the ENA metadata and files in a repeatable way.

The current production model is:

```text
one assembly + partial release
  └── one ENA Project/Study
        └── one ENA ANALYSIS per BAM/CRAM alignment file
```

Before project and analysis XML is generated, source run metadata is resolved
from the ENA Portal API. This supplies the authoritative BioSample, sample,
experiment, study, scientific-name, and taxon relationships. The resolved
file-level metadata is published under `ena_submission/metadata/` for audit.
Accession values supplied in the manifest are treated as expectations: valid
sample and experiment accession disagreements fail the run, while non-accession
labels are retained only as manifest provenance. XML is generated from the
resolved ENA values.

Aggregates use explicit source relationships: all represented source runs are
retained, and unique BioSamples are included when they can be represented
accurately. The workflow never creates a synthetic BioSample for a computational
aggregate. See [ENA_METADATA_SPEC.md](ENA_METADATA_SPEC.md) for the contract.

ENA does not accept multiple BAM/CRAM files in one `ANALYSIS`. The workflow
therefore expands the annotation-level manifest into one analysis per file,
while retaining the shared project metadata and source run/sample links.

## Why the workflow is split into these stages

1. **Build the manifest** from the RNA-seq annotation directory. This records
   which alignments exist and the metadata needed to describe them.
2. **Prepare the reference** once. A FASTA index and, when required, reference
   supplements are shared by all alignment tasks.
3. **Validate and reheader alignments** so BAM/CRAM sequence names and lengths
   agree with the INSDC reference used for the submission.
4. **Convert and index when requested.** BAM inputs can be converted to CRAM;
   the submitted checksums and indexes are calculated from the final files.
5. **Generate ENA XML and upload files.** Project and analysis XML are created
   from the manifests, while alignment files are uploaded to the Webin dropbox.
6. **Poll for accessions and publish run information.** The final output links
   submitted analyses to ENA accessions and retains pipeline reports for audit.

This separation keeps metadata decisions independent from file preparation and
makes the expensive alignment work reusable across retries. TEST submissions
exercise the same metadata and upload path against ENA's test service before a
production submission is attempted.

For a complete first-time walkthrough, including TEST/PROD setup, credential
handling, manifest validation, output inspection, and troubleshooting, see
[WALKTHROUGH.md](WALKTHROUGH.md).

## Quick start

The checked-in example is an annotation-level manifest. From the repository
root, use it for a wiring-only test. This runs the pipeline in stub mode and
does not contact ENA:

```bash
nextflow run pipelines/ena_submit/main.nf \
  -c pipelines/ena_submit/nextflow.config \
  -profile test \
  -stub-run \
  --manifest pipelines/ena_submit/examples/manifest.tsv \
  --webin_user test-user \
  --reheader_bams false \
  --outdir /tmp/ena-submit-stub
```

For a real TEST or PROD submission, first store `ENA_WEBIN_PASSWORD` in the
Nextflow secrets store, then use a generated manifest:

```bash
read -rsp 'ENA Webin password: ' ENA_WEBIN_PASSWORD
printf '\n'
nextflow secrets set ENA_WEBIN_PASSWORD "$ENA_WEBIN_PASSWORD"
unset ENA_WEBIN_PASSWORD

nextflow run pipelines/ena_submit/main.nf \
  -c pipelines/ena_submit/nextflow.config \
  --manifest /path/to/ena_manifest/manifest.tsv \
  --mode test \
  --webin_user "$WEBIN_USER" \
  --outdir /path/to/ena_submit_results
```

Use `--mode prod` only after the TEST submission has been inspected and the
resulting metadata has been accepted. Never put a password on the command
line; the full credential setup is in [WALKTHROUGH.md](WALKTHROUGH.md).

### Add alignment analyses to an existing study

To submit new alignment analyses to a study that already exists and is owned by
the Webin account, put that study accession or alias in `study` and run with
`--create_project false`. The resolver will attach the authoritative existing
sample and run accessions from ENA to each analysis; it will not create or
duplicate samples.

Registering genuinely new ENA sample objects is a separate sample submission.
It requires a checklist-aware sample manifest containing taxon and sample
attributes, and should be completed before the associated analysis submission.

### Repair an existing analysis sample link

If an analysis was accepted without a `SAMPLE_REF`, retrieve its current XML
and prepare a Webin V1 `MODIFY` request. `SAMPLE_REF` must contain the ENA
Sample accession (`ERS`, `SRS`, or `DRS`), not the linked BioSample accession
(`SAMEA`, `SAMN`, or `SAMD`). The helper resolves the source `RUN_REF` through
the ENA Portal API, preserves the existing `FILE` declarations, and records
the BioSample accession as analysis-level provenance.

```bash
modify_dir=/tmp/ena-modify_ERZ30663408
mkdir -p "$modify_dir"

# This endpoint requires the Webin account that owns the analysis.
# `curl` prompts for the password without placing it in the process arguments.
curl --user "$WEBIN_USER" \
  "https://www.ebi.ac.uk/ena/submit/drop-box/analyses/ERZ30663408" \
  -o "$modify_dir/analysis.xml"

python3 pipelines/ena_submit/bin/generate_analysis_modify_xml.py \
  --analysis-xml "$modify_dir/analysis.xml" \
  --outdir "$modify_dir"
```

Inspect `analysis_modify.xml` before submitting. Validate the payload first:

```bash
curl --user "$WEBIN_USER" \
  -F "SUBMISSION=@$modify_dir/submission_validate.xml" \
  -F "ANALYSIS=@$modify_dir/analysis_modify.xml" \
  https://www.ebi.ac.uk/ena/submit/drop-box/submit/ \
  -o "$modify_dir/validation_receipt.xml"
```

Only after validation succeeds, submit the update with the `MODIFY` action:

```bash
curl --user "$WEBIN_USER" \
  -F "SUBMISSION=@$modify_dir/submission_modify.xml" \
  -F "ANALYSIS=@$modify_dir/analysis_modify.xml" \
  https://www.ebi.ac.uk/ena/submit/drop-box/submit/ \
  -o "$modify_dir/modify_receipt.xml"
```

The modify helper intentionally does not change analysis files. ENA treats a
file change as a new analysis submission; this repair is only for metadata
relationships. After ENA processes the request, confirm that the analysis XML
contains `SAMPLE_REF accession="ERS..."` and retains its original file path
and checksum.

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

- `manifest.tsv`: one row for the annotation-level release/project metadata.
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

Runs without a matching alignment are recorded in `missing_files.tsv` and reported as a loud warning; manifest creation continues with the alignments that are available. Use `--fail-on-missing-files` when a complete alignment set is required. The older `--allow-missing-files` option remains accepted for compatibility.

If `samtools` is available, the builder inspects BAM/CRAM headers for basic alignment provenance such as sort order, `@PG` programs, and STAR version. Use `--no-bam-header` to skip this.

## Step 2: Submit To ENA

### BAM reference-name preparation

By default, BAM files are reheadered before checksumming and upload so their
`@SQ SN:` values match the INSDC reference FASTA. CRAM files are validated
against the same reference but are not rewritten. The manifest builder
automatically finds the matching genomic FASTA and NCBI assembly report beside
the assembly directory. The workflow creates the FASTA `.fai` index once in a
reference-preparation task and shares it with the reheadering tasks.

Use these options only to override automatic manifest discovery:

```bash
--reference_fasta /path/to/insdc-reference.fa \
--reference_assembly_report /path/to/GCA_xxx_assembly_report.txt
```

The process validates every alignment reference name and sequence length. For
BAMs it preserves `@SQ` order and other header records, writes a new BAM, and
rebuilds its index. Set `--reheader_bams false` only when all input alignments
have already been validated against the exact INSDC reference.

Before submitting, store the Webin password in Nextflow's local secrets store. This
keeps it out of command-line parameters, `.command.sh` files, and task arguments:

```bash
read -rsp 'ENA Webin password: ' ENA_WEBIN_PASSWORD
printf '\n'
nextflow secrets set ENA_WEBIN_PASSWORD "$ENA_WEBIN_PASSWORD"
unset ENA_WEBIN_PASSWORD
```

The secret name must be exactly `ENA_WEBIN_PASSWORD`. The username is passed as
`--webin_user`; do not pass a password using `--webin_password`.

Real TEST and PROD runs require `ENA_WEBIN_PASSWORD` in the Nextflow secrets
store. Nextflow also checks declared secrets while constructing stub tasks, so
the wiring-only command above needs an isolated placeholder secret when run
locally or in CI:

```bash
export NXF_HOME="$(mktemp -d)"
nextflow secrets set ENA_WEBIN_PASSWORD stub-password
```

Use the same temporary `NXF_HOME` for the subsequent `nextflow run` command;
never use a production password for a stub run.

Run the Nextflow workflow with the generated annotation manifest:

```bash
nextflow run pipelines/ena_submit/main.nf \
  -c pipelines/ena_submit/nextflow.config \
  --manifest /path/to/ena_manifest/GCA_011064425.1/manifest.tsv \
  --mode test \
  --webin_user "$WEBIN_USER" \
  --umbrella_study PRJEB000000 \
  --umbrella_alias ensembl_alignment_evidence_umbrella \
  --umbrella_title "Ensembl alignment evidence umbrella study" \
  --umbrella_description "Umbrella study grouping Ensembl alignment evidence projects." \
  --release_mode none \
  --remote_dir GCA_011064425.1-Ensembl-2025-12 \
  --upload_parallelism 2 \
  --outdir /path/to/ena_submit_results/GCA_011064425.1
```

Use `--mode prod` for production endpoints after Webin test validation.

When adding a child project to an existing ENA umbrella study, provide the
umbrella accession, exact registered alias, title, and description. After the
child project and analyses have succeeded, the workflow submits a terminal
project-level `MODIFY` request containing the verified child project accession:

```bash
--umbrella_study PRJEB000000 \\
--umbrella_alias ensembl_alignment_evidence_umbrella \\
--umbrella_title "Ensembl alignment evidence umbrella study" \\
--umbrella_description "Umbrella study grouping Ensembl alignment evidence projects."
```

Release is opt-in. Use `--release_mode immediate` in production to release the
umbrella and child project, or `--release_mode date --release_date YYYY-MM-DD`
to set a common release date. Release actions are rejected in TEST mode.

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
  --outdir /path/to/ena_submit_results/GCA_011064425.1
```

The reference index must be at `/path/to/reference.fa.fai`. CRAM files and their
CRAI indexes are uploaded; only the CRAM files are listed in the ENA analysis XML.
Existing non-BAM inputs are passed through unchanged.

## Manifest Schema

`manifest.tsv` has one row per annotation/release:

- `files_tsv`: Path to the file-level manifest.
- `study`: Existing child study accession/alias. It must be blank when
  `--create_project true`; use `--create_project false` to add analyses to an
  existing study.
- `project_alias`: Child project alias, usually `prj_<assembly_accession>_Ensembl_<YYYY_MM>`.
- `umbrella_study`: Existing umbrella project accession. Enables the terminal
  project-level `MODIFY` link submission.
- `umbrella_alias`: Exact alias used when the umbrella project was registered.
- `umbrella_title`: Existing umbrella title; must match the registered object.
- `umbrella_description`: Existing umbrella description; must match the registered object.
- `analysis_alias`: Stable ENA analysis alias.
- `title`: Human-readable analysis title.
- `description`: Human-readable description.
- `assembly_accession`: INSDC assembly accession, e.g. `GCA_052040795.1`.
- `reference_fasta`: Absolute path to the discovered INSDC genomic FASTA.
- `assembly_report`: Absolute path to the discovered NCBI assembly report.
- `reference_supplement`: Absolute path to an automatically generated FASTA containing report-listed non-nuclear sequences absent from the genomic FASTA.
- `last_geneset_update`: Genome metadata value, e.g. `2025-12`.
- `partial_release_label`: Derived release label.
- `species`: Production species name.
- `taxon_id`: Optional taxon ID.
- `ref_seqs`: Automatically populated with supplemental INSDC sequence accessions such as `NC_012374.1` when the assembly report identifies them.
- `analysis_links`: Optional `Label|URL; Label2|URL2`.
- `analysis_attributes`: Optional `key=value; key2=value2`. The builder writes `attr_*` keys.
- `analysis_type`: `REFERENCE_ALIGNMENT`.
- `omit_run_refs_in_test`: `false` by default, so generated TEST analyses retain their source `RUN_REF` links. Set it to `true` when testing with runs that are absent from ENA TEST.

The run-manifest builder writes this nested `manifest.tsv` plus one `files/*.tsv`
per analysis. This is the format consumed by the workflow; the older flat
run-list layout is no longer sufficient.

`files.tsv` has one row per alignment file:

- `file_path`: Local BAM/CRAM path.
- `file_type`: `bam` or `cram`.
- `remote_name`: Filename to use in the Webin drop-box.
- `run_accession`: Source ENA/SRA/DRA run accession.
- `run_accessions`: Comma-separated source runs represented by one aggregate file.
- `sample_accession`: Source sample accession from the RNA-seq CSV.
- `sample_accessions`: Optional source samples for an aggregate; these are
  re-resolved from ENA before submission.
- `experiment_accession`: Optional experiment accession.
- `platform`: Sequencing platform from the RNA-seq CSV.
- `source_fastq_count`: Number of source FASTQ records for the run.
- `source_fastq_urls`: Comma-separated source FASTQ URLs.
- `source_fastq_md5s`: Comma-separated source FASTQ MD5s.
- `bam_sort_order`: Header-derived sort order when available.
- `bam_pg_programs`: Header-derived `@PG` program names when available.
- `alignment_software`: Header-derived aligner name when available.
- `alignment_software_version`: Header-derived aligner version when available.

The resolved metadata also retains secondary sample/study accessions, library
strategy/source, instrument platform/model, and a per-run JSON provenance record.
For registered ENA samples, `sample_accessions` contains the ENA Sample
accession (`ERS`/`SRS`/`DRS`) used in `SAMPLE_REF`; `biosample_accessions`
contains the linked BioSample accession (`SAMEA`/`SAMN`/`SAMD`).

Only values matching ENA sample accession formats are emitted as `SAMPLE_REF`
elements. Tissue labels or other descriptive values in the source CSV are ignored.
Aggregated alignments do not receive synthetic BioSamples: the resolver emits
all represented source samples when ENA provides them and retains the complete
run list and resolved metadata as provenance.

## Parameters

| Parameter | Required | Description |
| --- | --- | --- |
| `--manifest` | yes | Annotation-level manifest produced by the builder or supplied directly |
| `--webin_user` | yes | ENA Webin username; the password is read from Nextflow secrets |
| `--mode` | no | `test` by default; use `prod` for the production endpoint |
| `--outdir` | no | Published output directory, default `./results` |
| `--create_project` | no | Create a child project (`true`) or reuse the study in each manifest row (`false`). |
| `--reheader_bams` | no | Validate BAM/CRAM inputs against the INSDC reference and reheader BAMs, default `true` |
| `--convert_to_cram` | no | Convert BAM inputs to CRAM before upload, default `false` |
| `--reference_fasta` | conditional | Required for CRAM conversion unless discovered from the manifest |
| `--upload_parallelism` | no | Maximum number of concurrent uploads |

The complete parameter list and types are defined in
[`nextflow_schema.json`](nextflow_schema.json).

## Outputs

- `${outdir}/ena_submission/projects/`: generated child project XML.
- `${outdir}/ena_submission/metadata/`: ENA-resolved source relationships and summaries.
- `${outdir}/ena_submission/xml/`: generated analysis XML.
- `${outdir}/ena_submission/accessions.tsv`: Webin polling results.
- `${outdir}/pipeline_info/`: execution timeline, report, trace, and software versions.
- Nextflow work directory: upload logs and task details.

## Current Production Decisions

Locked in:

- One ENA `ANALYSIS` per alignment file within an annotation/assembly partial release.
- Include all source `RUN_REF` entries in production.
- Support BAM and CRAM input files.
- Derive partial release identity from assembly accession plus `last_geneset_update`.
- Resolve RNA-seq source relationships from the ENA Portal API, retaining the
  local RNA-seq metadata as an auditable input.
- Support explicit creation of a child project or reuse of an existing study
  with `--create_project`.
- Add analyses to an existing study by setting `--create_project false` and
  supplying that study in each manifest row. Existing source samples are
  referenced by accession; they are not copied or recreated.
- Support project-level umbrella updates using verified child accessions.
- Accept registry/core metadata later through CLI-compatible fields.

Pending:

- Whether production submissions should use BAMs as-is or CRAMs.
- Confirm service-level acceptance of cross-archive `SAMD...`/`SAMN...`
  identifiers as direct `SAMPLE_REF` values with representative Webin TEST
  submissions.
- Registering new sample objects is a separate ENA sample submission and is
  not inferred from alignment manifests. It requires a sample metadata/checklist
  input and must not create synthetic samples for aggregate alignments.
- Which registry fields will replace CLI arguments for assembly/annotation metadata.
