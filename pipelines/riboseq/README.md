# RiboSeq Pipeline

Ribosome profiling pipeline: organism setup, data acquisition/collapsing, QC, STAR alignment (genome+transcriptome), RiboMetric + RiboWaltz, BEDGraph/BigWig, optional TrackHub, and optional TranslonScorer scoring.

## Quick Start

```bash
nextflow run main.nf -stub --outdir results
```

## Reusing an existing run context

For a rerun, generate a small config from the reference bundle, sample sheet,
and existing collapsed-read directory. This validates the paths and records
the resolved inputs in `pipeline_info/run_context.json`; processing choices
remain explicit command-line options.

If the reference directory contains the generated `params.config`, it is read
automatically. Otherwise the script discovers the reference files below the
directory. A generic species parameter file containing only settings such as
`base` and `genome_version` is not a reference bundle and will be rejected;
point `--reference-config` at the generated config or use the directory that
contains the actual organism-setup outputs.

```bash
python3 pipelines/riboseq/bin/make_riboseq_run_config.py \
  --reference-dir "$REF" \
  --collapsed-read-path "$OLD_COLLAPSED" \
  --sample-sheet "$MANIFEST" \
  --outdir "$OUT" \
  --output "$OUT/pipeline_info/run_context.config"

export NXF_VER=26.04.6
nextflow run pipelines/riboseq/main.nf \
  -profile slurm \
  -resume \
  -work-dir "$WORK" \
  -c "$OUT/pipeline_info/run_context.config" \
  --fetch true \
  --force_fetch true \
  --rpf_extraction_method getrpf \
  --getrpf_command extract \
  --run_matrix_mode false
```

Use the same work directory for retries so Nextflow can reuse completed tasks.
The input-resolution table is published under
`$OUT/pipeline_info/input_resolution/`. Published sample outputs are not used
as pipeline inputs: `-resume` is the mechanism for task reuse, while
`storeDir` should be reserved for stable reference/cache data.

## Structure

```
.
├── main.nf                       # Subworkflow example (2 tools)
├── workflows/
│   ├── simple_workflow.nf       # Simplest: 1 process only
│   └── subworkflow_example.nf   # Same as main.nf (for reference)
├── subworkflows/                 # 2 reusable subworkflow examples
├── modules/                      # Pipeline processes (including TranslonScorer)
├── advanced_entrypoints/         # Advanced: Entry point system
└── docs/                         # All documentation
```

## Examples

### Simple Workflow (workflows/simple_workflow.nf)
A single process workflow - demonstrates the basics.

```bash
nextflow run workflows/simple_workflow.nf -stub --outdir results
```

**Output**: standard riboseq artifacts (see below)

### Main Workflow (main.nf)
Two subworkflows chained together (also available as `workflows/subworkflow_example.nf`).

```bash
nextflow run main.nf -stub --outdir results
```

**Output**: alignment, analysis, bedgraphs/bigwigs, optional trackhub and translonscorer outputs

### Advanced Entry Points (advanced_entrypoints/)
Dynamic entry point system with automatic workflow resumption.

```bash
cd advanced_entrypoints
nextflow run main.nf --help
nextflow run main.nf -stub --outdir results
```

**Output**: multiple entrypoints for development

## Where to Start

**New to Nextflow?** Start with `workflows/simple_workflow.nf` to understand the basics.

**Familiar with workflows?** Explore `main.nf` to see subworkflow patterns.

**Building complex pipelines?** Check out `advanced_entrypoints/` for the entry point system.

## Documentation

All documentation is in **[docs/](docs/)**:
- [docs/QUICK_START.md](docs/QUICK_START.md) - Complete guide
- [docs/PATTERNS.md](docs/PATTERNS.md) - Common patterns
- [docs/INDEX.md](docs/INDEX.md) - Full index

## Parameters

- `--outdir` - Output directory (default: `'results'`)

### TranslonScorer (optional)

Enable per‑sample scoring after post‑processing for QC‑passing samples:

```bash
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet <samples.csv> \
  --star_index <STAR_index> \
  --gtf <anno.gtf> \
  --fasta <genome.fa> \
  --chrom_sizes_file <chrom.sizes> \
  --run_translonscorer true \
  --translonscorer_bam_type unique_no_junction \
  -profile local -stub
```

Outputs per sample: `<outdir>/translonscorer/<sample>/<sample>_orfs_scored.csv`, `<sample>_report.html`.

### ChOROS sequence-bias correction pilot

ChOROS can be run once per sample after the `good` QC gate. It consumes the
STAR transcriptome BAM and RiboMetric A-site offsets, preserves collapsed-read
abundance through the `ZW` BAM tag, and publishes corrected transcript/codon
counts without replacing the baseline BigWigs.

Build and publish a pinned runtime from
`containers/choros/Dockerfile`, then enable the stage with:

```bash
nextflow run pipelines/riboseq/main.nf \
  --sample_sheet <samples.csv> \
  --star_index <STAR_index> \
  --gtf <annotation.gtf> \
  --fasta <genome.fa> \
  --transcriptome_fasta <transcripts.fa> \
  --ribometric_annotation <ribometric.tsv> \
  --chrom_sizes_file <chrom.sizes> \
  --run_choros true \
  --choros_container <registry>/choros:23d3e424 \
  -profile slurm
```

Organism setup supplies `transcriptome_fasta` and `ribometric_annotation`
automatically. ChOROS requires `ribometric_offset_target = 'a_site'`. Outputs
are written under `<outdir>/choros/<sample>/`; the BAM preparation metrics
report how many query groups survived the conservative single-transcript rule.
