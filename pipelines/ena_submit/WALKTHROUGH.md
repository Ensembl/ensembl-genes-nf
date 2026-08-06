# ENA alignment submission walkthrough

This walkthrough takes a new user from an Ensembl genebuild RNA-seq directory to an ENA submission.

The workflow creates:

\`\`\`text
one Project/Study per assembly + partial release
  └── one ANALYSIS per BAM or CRAM alignment file
\`\`\`

ENA does not accept multiple BAM/CRAM files in one ANALYSIS. TEST is safe for validation, but TEST objects are temporary and removed by ENA.

## 1. Prerequisites

You need:

- this repository on the \`feature/ena_submission\` branch, or a later branch containing the ENA pipeline;
- Nextflow and Java;
- Singularity/Apptainer;
- access to the RNA-seq directory and alignment files;
- an ENA Webin account;
- the cluster Slurm profile.

Check the tools:

\`\`\`bash
nextflow -version
java -version
singularity --version
\`\`\`

Pull the implementation:

\`\`\`bash
git switch feature/ena_submission
git pull --ff-only
git log -1 --oneline
\`\`\`

## 2. Store the Webin password safely

Do not pass the password as \`--webin_password\`. That exposes it in shell history, Nextflow commands, and task files.

\`\`\`bash
read -rsp 'ENA Webin password: ' P
printf '\\n'
nextflow secrets set ENA_WEBIN_PASSWORD "$P"
unset P
nextflow secrets list
\`\`\`

The secret name must be exactly \`ENA_WEBIN_PASSWORD\`. Pass only the username normally:

\`\`\`bash
export WEBIN_USER='Webin-XXXXX'
\`\`\`

If a password has appeared in a transcript or log, rotate it through ENA.

## 3. Check the input layout

The builder expects:

\`\`\`text
<species>/<assembly_accession>/rnaseq/
├── <species>.csv
└── output/
    ├── <RUN>_Aligned.sortedByCoord.out.bam
    └── ...
\`\`\`

It assumes a headerless tab-delimited CSV where column 1 is the source sample field, column 2 is the run accession, column 9 is the platform, column 11 contains FASTQ URLs, and column 12 contains FASTQ MD5 values.

Set and inspect the rat input:

\`\`\`bash
RNASEQ=/hps/nobackup/flicek/ensembl/genebuild/ereboperezsilva/annot_vert/rattus_rattus/GCA_011064425.1/rnaseq

find "$RNASEQ/output" -maxdepth 1 -type f -name '*.bam' | wc -l
find "$RNASEQ/output" -maxdepth 1 -type f \( -name '*.bam' -o -name '*.cram' \) | head
\`\`\`

If files are nested, inaccessible, or named differently, fix discovery before treating them as missing. Missing files are reported as warnings by default; use \`--fail-on-missing-files\` when a complete set is required.

## 4. Build a TEST manifest

Use a new alias namespace for each TEST run because ENA removes TEST objects.

\`\`\`bash
ASSEMBLY=GCA_011064425.1
RELEASE=2025-12
MANIFEST_ROOT=/hps/nobackup/flicek/ensembl/genebuild/jackt/ena-submission/ena_manifest

python3 pipelines/ena_submit/bin/build_manifest_from_rnaseq_annotation.py \
  --rnaseq-dir "$RNASEQ" \
  --assembly-accession "$ASSEMBLY" \
  --last-geneset-update "$RELEASE" \
  --species rattus_rattus \
  --file-format bam \
  --allow-missing-files \
  --project-alias prj_GCA_011064425.1-Ensembl-2025-12-test-20260728 \
  --analysis-alias rnaseq_alignment_evidence_GCA_011064425.1_Ensembl_2025-12_test_20260728 \
  --outdir "$MANIFEST_ROOT/GCA_011064425.1_test_20260728"
\`\`\`

Inspect the generated files:

\`\`\`bash
OUT_MANIFEST="$MANIFEST_ROOT/GCA_011064425.1_test_20260728"

column -t -s $'\\t' "$OUT_MANIFEST/summary.tsv"
column -t -s $'\\t' "$OUT_MANIFEST/missing_files.tsv" | head -30
column -t -s $'\\t' "$OUT_MANIFEST/files.tsv" | head -10
\`\`\`

The rat example previously contained 99 metadata runs and 52 matching alignment files. The missing runs are recorded explicitly and should be understood before production.

## 5. Run the TEST submission

\`\`\`bash
RESULTS=/hps/nobackup/flicek/ensembl/genebuild/jackt/ena-submission/ena_results/GCA_011064425.1_test_20260728

nextflow run pipelines/ena_submit/main.nf \
  -c pipelines/ena_submit/nextflow.config \
  -profile singularity,slurm \
  --manifest "$OUT_MANIFEST/manifest.tsv" \
  --mode test \
  --webin_user "$WEBIN_USER" \
  --remote_dir GCA_011064425.1-Ensembl-2025-12-test-20260728 \
  --upload_parallelism 2 \
  --outdir "$RESULTS" \
  -resume
\`\`\`

Never add \`--webin_password\`.

The stages are:

\`\`\`text
project XML → submit/poll project
  → expand manifest into one row per alignment
  → MD5 → FTP upload
  → one analysis XML per alignment
  → submit/poll analyses
\`\`\`

For the rat run, success should show 52 tasks for MD5, FTP, XML generation, and analysis submission. The analysis polling process may show \`1 of 1\` because one task polls all 52 queue receipts.

## 6. Inspect the result locally

\`\`\`bash
find "$RESULTS/ena_submission" -type f | sort
column -t -s $'\\t' "$RESULTS/ena_submission/accessions.tsv"
awk -F '\\t' 'NR > 1 { n[$4]++ } END { for (x in n) print x, n[x] }' \
  "$RESULTS/ena_submission/accessions.tsv"
\`\`\`

Inspect one analysis XML:

\`\`\`bash
ANALYSIS_XML=$(find "$RESULTS/ena_submission/xml" -name analysis.xml | head -1)
xmllint --format "$ANALYSIS_XML" | less
\`\`\`

Every analysis should contain exactly one alignment file:

\`\`\`bash
for f in "$RESULTS"/ena_submission/xml/*/analysis.xml; do
  printf '%s: ' "$f"
  grep -c '<FILE ' "$f"
done
\`\`\`

In TEST mode, RUN_REF is omitted by default because source runs may only exist in PROD. Tissue labels and other non-accession sample values are ignored.

## 7. View the TEST submission in ENA

While the records still exist, use:

- [TEST Webin Portal](https://wwwdev.ebi.ac.uk/ena/submit/webin/login)
- [TEST Webin Reports](https://wwwdev.ebi.ac.uk/ena/submit/report)

Search using accessions from \`accessions.tsv\`. The portal can show metadata, submitted files, archival status, and processing status. TEST records are discarded, so keep the local XML, queue receipts, and accession table.

## 8. Move from TEST to PROD

Before PROD, confirm:

- the assembly and release metadata are correct;
- every analysis XML has exactly one FILE;
- remote filenames are unique;
- the partial-data decision is intentional;
- sample values are real ENA accessions or omitted;
- the umbrella study is correct;
- production Webin credentials are available.

Create a fresh production manifest with production aliases:

\`\`\`bash
python3 pipelines/ena_submit/bin/build_manifest_from_rnaseq_annotation.py \
  --rnaseq-dir "$RNASEQ" \
  --assembly-accession "$ASSEMBLY" \
  --last-geneset-update "$RELEASE" \
  --species rattus_rattus \
  --file-format bam \
  --allow-missing-files \
  --project-alias prj_GCA_011064425.1-Ensembl-2025-12 \
  --analysis-alias rnaseq_alignment_evidence_GCA_011064425.1_Ensembl_2025-12 \
  --outdir "$MANIFEST_ROOT/GCA_011064425.1_prod"
\`\`\`

Launch PROD with a new output directory and without reusing TEST state:

\`\`\`bash
PROD_RESULTS=/hps/nobackup/flicek/ensembl/genebuild/jackt/ena-submission/ena_results/GCA_011064425.1_prod

nextflow run pipelines/ena_submit/main.nf \
  -c pipelines/ena_submit/nextflow.config \
  -profile singularity,slurm \
  --manifest "$MANIFEST_ROOT/GCA_011064425.1_prod/manifest.tsv" \
  --mode prod \
  --webin_user "$WEBIN_USER" \
  --remote_dir GCA_011064425.1-Ensembl-2025-12 \
  --upload_parallelism 2 \
  --outdir "$PROD_RESULTS"
\`\`\`

PROD creates permanent ENA records. Keep local copies of all uploaded files.

## 9. Optional BAM-to-CRAM conversion

CRAM conversion is opt-in and requires a matching reference FASTA and index:

\`\`\`bash
--convert_to_cram true \
--reference_fasta /path/to/reference.fa
\`\`\`

The index must be \`/path/to/reference.fa.fai\`. The workflow converts before MD5 calculation, creates CRAI indexes, uploads CRAM plus CRAI, and lists only the CRAM in the analysis XML.

## 10. Troubleshooting

### \`530 Login incorrect\`

Reset the secret without putting the password in the command line and test interactively:

\`\`\`bash
read -rsp 'ENA Webin password: ' P
printf '\\n'
nextflow secrets set ENA_WEBIN_PASSWORD "$P"
unset P
lftp -u "$WEBIN_USER" webin2.ebi.ac.uk
\`\`\`

### \`Invalid group of files: 52 "bam" files\`

The old annotation-level XML or a stale cached expansion was used. Pull the current branch and confirm analysis aliases end with run accessions such as \`_DRR503133\`.

### \`Failed to find referenced study\`

A TEST project from an earlier day was discarded. Use a new project alias and a new TEST output directory.

### \`Failed to find referenced sample, accession "ileum"\`

The source contains a tissue label rather than an ENA sample accession. Pull the latest implementation, which ignores non-accession sample values.

### Slurm says \`Neither --mem nor --mem-per-cpu specified\`

Use \`-profile singularity,slurm\` and pull the branch containing shared resource configuration.

### \`47 runs had no alignment file\`

Check \`missing_files.tsv\`, directory permissions, and the naming pattern. Missing files are warnings by default; use \`--fail-on-missing-files\` when the partial set is not intentional.

## 11. What to retain

Keep these items together:

- the exact Git commit;
- the manifest builder command;
- manifest.tsv, files.tsv, missing_files.tsv, and summary.tsv;
- the Nextflow command with the password removed;
- project and analysis XML files;
- accessions.tsv;
- the final Nextflow summary and ENA receipts.

Never retain passwords in presentations, shell history, task logs, or shared output directories.
