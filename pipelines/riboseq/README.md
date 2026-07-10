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
