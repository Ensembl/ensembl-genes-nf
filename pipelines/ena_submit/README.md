# ENA Submission Pipeline

This pipeline uploads analysis files (e.g. BAM/CRAM) to ENA Webin drop-box and submits an ANALYSIS linking them to public runs.
See also: `pipelines/ena_submit/SUBMISSION_STRATEGY.md` for project/alias and idempotency conventions.

## Inputs

Provide a tab-separated `manifest.tsv` with the following columns (header required):

- file_path: Absolute or relative path to the file to submit (BAM/CRAM).
- file_type: One of `bam`, `cram`.
- study: ENA Study accession or alias (e.g. PRJEB12345 or your alias).
- analysis_alias: Unique alias for the analysis (string). If empty, the pipeline generates one.
- title: Human-readable title.
- description: Short description.
- run_accessions: Comma-separated ENA run accessions that the file derives from (e.g. ERR123,ERR456).
- run_list_path: Optional path to a text/TSV file containing run accessions (commas and whitespace allowed). Merged with `run_accessions` and de-duplicated.
- experiment_accessions: Optional comma-separated ENA experiment accessions (if you prefer linking experiments).
- assembly_accession: ENA/INSDC assembly accession the reads were aligned to (e.g. GCA_000001405.28 or GCF_...).
- sample_accession: Optional sample accession to associate (e.g. ERS1234567). Leave blank for cross-sample merges.
- ref_seqs: Optional comma-separated list of reference sequence accessions (for <SEQUENCE> entries) when relevant.
- remote_name: Optional alternate remote filename (defaults to the source basename).
- analysis_type: Optional; one of `READ_ALIGNMENT` or `REFERENCE_ALIGNMENT` (default). The pipeline coerces `READ_ALIGNMENT` to `REFERENCE_ALIGNMENT` to satisfy current ENA schema validation.
- analysis_links: Optional; `Label|URL; Label2|URL2`.
- analysis_attributes: Optional; `key=value; key2=value2`. You can also supply `attr_*` columns; e.g., `attr_pipeline=ensembl-genes-nf` becomes a TAG/VALUE pair.
- omit_run_refs_in_test: Optional; `true|false` to omit or emit RUN_REF when `--mode test` (default true/omit).

Example: see `pipelines/ena_submit/examples/manifest.tsv`.

## Credentials

Set Webin credentials via environment variables before running:

- `WEBIN_USER` (e.g. `Webin-XXXXX`)
- `WEBIN_PASSWORD`

## Running

```
nextflow run pipelines/ena_submit/main.nf \
  --manifest pipelines/ena_submit/examples/manifest.tsv \
  --mode test \
  --release Ensembl_110 \
  --upload_protocol aspera \
  --upload_parallelism 2 \
  --ascp_limit 300M \
  --remote_dir myproject/subset1 \
  --submit_api v1 \
  --outdir results
```

Notes:
- Use `--mode prod` to switch to production endpoints.
- Use `--upload_protocol ftp` if Aspera is unavailable.
- Keep `--upload_parallelism` low to avoid overwhelming ENA (2–3 typical).
- In TEST mode, the pipeline omits RUN_REF by default so ENA TEST doesn’t 404 production run accessions. Set `omit_run_refs_in_test=false` in your manifest row to include RUN_REF in TEST.
- Projects are auto‑derived and registered per assembly+release; supply `--release`.

## Outputs

- `${outdir}/ena_submission/uploads/` – logs of uploads
- `${outdir}/ena_submission/xml/` – generated XML per analysis (per-ID subfolders)
- `${outdir}/ena_submission/receipts/` – ENA receipts (XML) per analysis
- `${outdir}/pipeline_info/` – Nextflow execution reports

## Submission model

- One ANALYSIS per file (ENA permits a single BAM/CRAM per READ/REFERENCE_ALIGNMENT).
- Files are uploaded first; the pipeline then submits metadata referencing the file with MD5.
- The pipeline can use Webin REST drop-box (v1) or Webin REST v2 (async).

## Test vs prod

- `--mode test` uses ENA test endpoints for submission. Files still upload to your Webin drop-box.
- `--mode prod` uses production endpoints.

## Safety

- Upload concurrency is capped with `--upload_parallelism`.
- Aspera uploads are bandwidth-limited with `--ascp_limit` (default 300M).

## Requirements

- For Aspera: `ascp` CLI available and Webin SSH access.
- For FTP: `lftp` installed.
- For submission: `curl` installed.
- Credentials exported: `WEBIN_USER`, `WEBIN_PASSWORD`.

## Utilities

- Build a manifest from a runs/tissue CSV (like `tissueall.csv`):

```
python3 pipelines/ena_submit/bin/build_manifest_from_runs.py \
  --input tissueall.csv \
  --study PRJEB999999 \
  --species-map species.tsv \
  --release Ensembl_110 \
  --group-by assembly,tissue \
  --file-dir /path/to/merged/files \
  --file-ext cram \
  --outdir /tmp/ena_manifest \
  --links-prefix https://example.org/runs
```

This writes `/tmp/ena_manifest/manifest.tsv` and per-group run lists under `/tmp/ena_manifest/runs/`. The manifest is compatible with this pipeline and supports per‑run mode by default. In TEST, `RUN_REF` are omitted unless overridden; in PROD they are included so your analyses link back to the submitter’s RUN and BioSample where available.
