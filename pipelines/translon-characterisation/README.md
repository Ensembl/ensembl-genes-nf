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

## Slurm execution

Keep the repository, Nextflow state, work directory, temporary files, container
cache, references, and results below a shared no-backup root.  The `slurm`
profile uses the Slurm executor and Singularity/Apptainer with automatic input
mounting.  Set `NXF_SINGULARITY_CACHEDIR` explicitly so the inherited shared
configuration does not place images elsewhere.

```bash
export TRANSLON_ROOT=/hps/nobackup/flicek/ensembl/genebuild/jackt/translon-characterisation
export NXF_HOME="$TRANSLON_ROOT/nextflow"
export NXF_SINGULARITY_CACHEDIR="$TRANSLON_ROOT/singularity-cache"
export NXF_TEMP="$TRANSLON_ROOT/tmp"
export TMPDIR="$TRANSLON_ROOT/tmp"

mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$NXF_TEMP" \
  "$TRANSLON_ROOT/work" "$TRANSLON_ROOT/results" "$TRANSLON_ROOT/logs"

cd "$TRANSLON_ROOT/repo/ensembl-genes-nf/pipelines/translon-characterisation"

nextflow run . \
  -profile slurm \
  -work-dir "$TRANSLON_ROOT/work" \
  --intervals "$TRANSLON_ROOT/inputs/annotated_orfs.bed" \
  --gencode_gff3 "$TRANSLON_ROOT/references/gencode.annotation.gff3" \
  --genome_fasta "$TRANSLON_ROOT/references/GRCh38.primary_assembly.genome.fa" \
  --translation_verdicts "$TRANSLON_ROOT/inputs/translonscorer.tsv" \
  --proteome_fasta "$TRANSLON_ROOT/references/gencode.pc_translations.fa" \
  --phylocsf_matched_null "$TRANSLON_ROOT/references/phylocsf_58mammals_null.json" \
  --outdir "$TRANSLON_ROOT/results" \
  -with-report "$TRANSLON_ROOT/logs/execution-report.html" \
  -with-trace "$TRANSLON_ROOT/logs/execution-trace.tsv" \
  -with-timeline "$TRANSLON_ROOT/logs/execution-timeline.html" \
  -resume
```
