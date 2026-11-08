# Translon characterisation

This pipeline characterises trusted translons; it does not rediscover ORFs or
recompute translation evidence.  The annotation unit is an `(interval ×
transcript)` instance.  Every axis returns `positive`, `null`,
`uninformative`, or `unattributable`; there is intentionally no aggregate
score.

## Current vertical slice

`SUBSTRATE_TOPOLOGY`, `TRANSLATION_JOIN`, typed-axis contracts, late set-level
clustering barrier, and `ADJUDICATE` are executable.  The initial Python
implementation establishes coordinate normalization, transcript-contingent
instances, soft retention of unhosted loci, and the claim wall.  Commodity
tool stages currently emit an honest `uninformative` state until their pinned
wrapper and reference inputs are provided; they must never manufacture a
negative result from absent data.

Run (from repository root):

```bash
nextflow run pipelines/translon-characterisation \
  --intervals intervals.tsv --gencode_gff3 gencode.gff3 \
  --translation_verdicts translonscorer.tsv --proteome_fasta proteins.fa \
  -profile local
```

Input intervals are headered TSV with `interval_id`, `chrom`, `start`, `end`,
`strand`, and optional `frame`; coordinates are 0-based, half-open.

## Optional evolutionary-history candidate evidence

The pipeline can attach a read-only per-species candidate record to each
attributed instance. Pass a JSONL or TSV file with one row per candidate and a
`species` plus either `instance_id` or `interval_id`. Optional fields include
`candidate_gene`, `candidate_transcript`, `orthology_status`, `synteny_status`,
`cds_status`, `leader_status`, `orf_status`, `translation_status`, and
`path_status`, and `orthology_source`. When source coordinates (`chrom`,
`start`, `end`, `strand`) are present, basic path validity is checked without
changing the candidate rows. Reference identifiers (`reference_gene`,
`reference_transcript`, `reference_cds`, and `reference_uorf`) are preserved
when present:

```bash
nextflow run pipelines/translon-characterisation \
  --intervals intervals.tsv --gencode_gff3 gencode.gff3 \
  --translation_verdicts translonscorer.tsv --proteome_fasta proteins.fa \
  --genome_fasta genome.fa --evolutionary_candidates candidates.jsonl \
  -profile local
```

This stage does not select an ortholog, calculate conservation, or change an
ORBL/BLS denominator. It preserves one-to-many candidates and emits an
`uninformative` typed axis until those candidates have passed a separate path
validation and scoring stage. A descriptive species summary is included, but
the candidate rows remain authoritative when a summary is ambiguous.

## Slurm execution

Keep the repository, work directory, temporary files, container cache, references,
and results below a shared no-backup root. The `slurm_apptainer` profile submits
child jobs, while `apptainer_local` runs inside an existing `salloc`; both use
automatic input mounting. The generated config pins the cache location outside
the work directory.

The setup helper creates this layout idempotently, checks Nextflow/Java and the
container runtime, and writes an explicit `hpc.config` fragment. It does not
download or copy large references, and it never submits a job:

```bash
python3 pipelines/translon-characterisation/bin/setup_hpc_environment.py \
  --root /hps/nobackup/flicek/ensembl/genebuild/jackt/translon-characterisation \
  --repo /hps/nobackup/flicek/ensembl/genebuild/jackt/ensembl-genes-nf \
  --executor local
```

If the inputs are already staged, add `--inputs-manifest inputs.tsv` to check
them before creating the layout. The manifest columns are `name`, `path`, and
optional `required` (`true` by default).

Then run a normal Nextflow command from your `salloc`, using absolute paths.
The generated config supplies the executor and Apptainer cache; input paths
remain explicit command-line arguments. With `--executor local`, tasks run
inside the resources already granted by `salloc`. Use `--executor slurm` when
you want Nextflow to submit child jobs instead:

```bash
nextflow run /hps/nobackup/.../ensembl-genes-nf/pipelines/translon-characterisation \
  -profile apptainer_local \
  -c /hps/nobackup/.../translon-characterisation/hpc.config \
  -work-dir /hps/nobackup/.../translon-characterisation/work \
  --intervals /hps/nobackup/.../inputs/annotated_orfs.bed \
  --gencode_gff3 /hps/nobackup/.../references/gencode.annotation.gff3 \
  --genome_fasta /hps/nobackup/.../references/GRCh38.primary_assembly.genome.fa \
  --translation_verdicts /hps/nobackup/.../inputs/translonscorer.tsv \
  --proteome_fasta /hps/nobackup/.../references/gencode.pc_translations.fa \
  --evolutionary_candidates /hps/nobackup/.../inputs/evolutionary_candidates.jsonl \
  --outdir /hps/nobackup/.../translon-characterisation/results \
  -with-report /hps/nobackup/.../translon-characterisation/logs/execution-report.html \
  -with-trace /hps/nobackup/.../translon-characterisation/logs/execution-trace.tsv \
  -with-timeline /hps/nobackup/.../translon-characterisation/logs/execution-timeline.html \
  -resume
```

Use `--check-only` first on a new cluster. The helper requires Nextflow 26.04+,
Java 17 or 21, and Apptainer/Singularity; it only warns when `sbatch` is not
available so the same check can be run on a preparation node.
