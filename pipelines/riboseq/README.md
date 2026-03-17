# RiboSeq Pipeline

Ribosome profiling pipeline: organism setup, data acquisition/collapsing, QC, STAR alignment (genome+transcriptome), RiboMetric + RiboWaltz, BEDGraph/BigWig, optional TrackHub, and optional TranslonScorer scoring.

## Quick Start

```bash
nextflow run main.nf -stub --outdir results
```

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
