# Validation Framework

This directory contains configuration files for systematic validation of different pipeline parameters.

## Quick Start

Run the pipeline with a specific config:
```bash
nextflow run main.nf \
  -params-file validation/configs/baseline.yaml \
  --sample_sheet samples.csv \
  --star_index /path/to/star \
  --gtf /path/to/annotation.gtf \
  --fasta /path/to/genome.fa \
  --chrom_sizes_file /path/to/chrom.sizes \
  --outdir results/baseline
```

## Configuration Files

### Baseline
- **baseline.yaml** - Default configuration for comparison baseline

### Adapter Detection & RPF Extraction

| Config | Method | Description |
|--------|--------|-------------|
| baseline.yaml | getrpf (extract) | Alignment-based RPF extraction (default) |
| adapter_getrpf_extract.yaml | getrpf extract | Explicit alignment-based extraction |
| adapter_getrpf_extractrpf.yaml | getrpf extract-rpf | Pattern matching + HMM segmentation |
| adapter_getrpf_decidetrim.yaml | getrpf decide-trim | Combined architecture + alignment decision |
| adapter_fastp_auto.yaml | fastp auto-detect | fastp built-in adapter detection |
| adapter_fastp_reflist.yaml | fastp + ref list | fastp with reference adapter list |
| rpf_traditional.yaml | traditional | FASTQC + fastp trimming |

### rRNA Filtering

| Config | Method | Reference | Description |
|--------|--------|-----------|-------------|
| baseline.yaml | None | - | No rRNA filtering (default) |
| rrna_filter_bowtie.yaml | Bowtie | GTF-extracted | Bowtie alignment vs Ensembl rRNAs |
| rrna_bowtie_ensembl.yaml | Bowtie | Ensembl GTF | Explicit Ensembl rRNA filtering |
| rrna_bowtie_silva.yaml | Bowtie | SILVA | Bowtie alignment vs SILVA database |
| rrna_ribodetector.yaml | RiboDetector | ML-based | Machine learning detection (no ref) |

### Alignment Parameters

#### Mismatches
| Config | Mismatches | Description |
|--------|------------|-------------|
| mismatch_0.yaml | 0 | Exact matches only |
| mismatch_1.yaml | 1 | 1 mismatch allowed |
| mismatch_2.yaml | 2 | 2 mismatches allowed |
| baseline.yaml | 3 | 3 mismatches (default) |

#### Alignment Mode & Soft-clipping
| Config | Mode | Clipping | Description |
|--------|------|----------|-------------|
| baseline.yaml | EndToEnd | None | No soft-clipping (default) |
| align_local.yaml | Local | Both ends | Full soft-clipping enabled |
| align_local_clip5.yaml | Local | 5' trim | 3bp 5' trimming + soft-clip |
| align_local_clip3.yaml | Local | 3' trim | 3bp 3' trimming + soft-clip |

### Multi-mapper Handling

| Config | Max Locations | Description |
|--------|---------------|-------------|
| multimap_1.yaml | 1 | Unique mappers only |
| multimap_5.yaml | 5 | Up to 5 locations |
| baseline.yaml | 10 | Up to 10 locations (default) |
| multimap_20.yaml | 20 | Up to 20 locations |

### Read Length Filtering

| Config | Range | Description |
|--------|-------|-------------|
| baseline.yaml | None | No filtering (default) |
| length_filter_strict.yaml | 28-32nt | Strict RPF range |
| length_filter_relaxed.yaml | 25-35nt | Relaxed RPF range |

### P-site Offset Calculation

#### RiboWaltz Methods
| Config | Extremity | Description |
|--------|-----------|-------------|
| baseline.yaml | auto | Automatic selection (default) |
| psite_5end.yaml | 5end | Calculate from 5' end |
| psite_3end.yaml | 3end | Calculate from 3' end |

#### RiboMetric Methods
| Config | Method | Description |
|--------|--------|-------------|
| baseline.yaml | tripsviz | TripsViz-style calculation (default) |
| offset_ribometric_changepoint.yaml | changepoint | Changepoint detection |
| offset_ribometric_ribowaltz.yaml | ribowaltz | RiboWaltz-style in RiboMetric |
| offset_external_ribowaltz.yaml | external | Pass RiboWaltz offsets to RiboMetric |
| offset_global_12.yaml | global | Fixed 12nt offset |
| offset_global_15.yaml | global | Fixed 15nt offset |

## Running Validation Suite

To run multiple configurations for comparison, use a simple loop:

```bash
# Define configs to test
CONFIGS=(
  # Baseline
  baseline

  # Adapter detection
  adapter_getrpf_extract
  adapter_getrpf_extractrpf
  adapter_getrpf_decidetrim
  adapter_fastp_auto

  # rRNA removal
  rrna_bowtie_ensembl
  rrna_bowtie_silva
  rrna_ribodetector

  # Alignment stringency
  mismatch_0
  mismatch_1
  mismatch_2
  align_local
  align_local_clip5
  align_local_clip3

  # Multi-mapper handling
  multimap_1
  multimap_5
  multimap_20

  # Read length filtering
  length_filter_strict
  length_filter_relaxed

  # P-site offset methods
  psite_5end
  psite_3end
  offset_ribometric_changepoint
  offset_ribometric_ribowaltz
  offset_external_ribowaltz
  offset_global_12
  offset_global_15
)

# Run each config
for config in "${CONFIGS[@]}"; do
  nextflow run main.nf \
    -params-file validation/configs/${config}.yaml \
    --sample_sheet samples.csv \
    --star_index /path/to/star \
    --gtf /path/to/annotation.gtf \
    --fasta /path/to/genome.fa \
    --chrom_sizes_file /path/to/chrom.sizes \
    --outdir results/${config} \
    -profile slurm \
    -resume
done
```

## Metrics to Compare

Key metrics for validation analysis:

### Read Processing
- Total reads (raw)
- Reads after collapsing
- Reads after rRNA filtering (if enabled)
- Aligned reads (unique vs multi-mapped)

### Alignment Quality
- Unique mapping rate
- Multi-mapping rate
- Unmapped rate
- Soft-clip length distribution (for Local mode)

### Frame Distribution
- Frame 0 percentage (should be ~65-70% for good data)
- Frame 1 and 2 percentages
- Distribution by read length

### Metagene Profiles
- Start codon enrichment
- Stop codon profile shape
- Periodicity signal strength

### Coverage
- CDS coverage depth
- 5'UTR/3'UTR ratio

### Offset Quality
- Offset by read length (RiboWaltz)
- Offset confidence scores
- Comparison between RiboWaltz and RiboMetric offsets

## Output Structure

Each config run creates:
```
results/{config_name}/
├── collapsed/           # Collapsed FASTA files
├── fastqc/              # Initial QC reports
├── rrna_filter/         # rRNA filtering logs (if enabled)
├── alignment/           # BAM files and stats
├── RiboMetric/          # RiboMetric QC reports
├── ribowaltz/           # P-site analysis
├── bedgraphs/           # Coverage tracks
├── bigwigs/             # BigWig files
└── pipeline_info/       # Nextflow reports
```

## Creating Custom Configs

Copy baseline.yaml and modify specific parameters:

```bash
cp validation/configs/baseline.yaml validation/configs/my_custom.yaml
# Edit my_custom.yaml with your parameters
```

Key parameter categories:
- RPF extraction: `rpf_extraction_method`, `getrpf_command`, `getrpf_max_reads`
- Adapter detection: `fastp_detect_adapter`, `fastp_adapter_sequence`, `fastp_adapter_fasta`
- rRNA filtering: `run_rrna_filter`, `rrna_filter_method`, `bowtie_rrna_mismatches`, `ribodetector_len`
- Alignment: `mismatches`, `max_multimappers`, `alignment_type`, `clip_5p_nbases`, `clip_3p_nbases`
- Length filtering: `filter_by_length`, `rpf_length_min`, `rpf_length_max`
- RiboMetric offset: `ribometric_offset_method`, `ribometric_offset_global`, `ribometric_use_ribowaltz_offsets`
- RiboWaltz: `ribowaltz_extremity`, `ribowaltz_flanking`

## Full Validation Matrix

For the complete validation matrix with all dimensions, metrics, and visualization requirements, see [VALIDATION_MATRIX.md](VALIDATION_MATRIX.md).
