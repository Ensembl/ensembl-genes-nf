# Ribo-Seq Validation Matrix

## Purpose

This document outlines a comprehensive validation framework to systematically evaluate Ribo-seq processing choices. The goal is to generate publication-ready evidence that justifies our final workflow configuration.

---

## Validation Dimensions

### 1. Adapter Detection & Removal

| Config ID | Method | Description | Key Params |
|-----------|--------|-------------|------------|
| `adapter_fastqc_tophit` | FastQC top hit | Use top overrepresented sequence from FastQC | - |
| `adapter_fastp_auto` | fastp auto-detect | Let fastp auto-detect adapters | - |
| `adapter_fastp_reflist` | fastp reference list | Trim using reference adapter list | `adapter_list` |
| `adapter_getrpf_extract` | getRPF extract | Alignment-based RPF extraction (recommended) | `sample_size`, `preserve_umi` |
| `adapter_getrpf_decidetrim` | getRPF decide-trim | Combined architecture + alignment decision | `max_reads` |

**getRPF Commands Available:**
- `extract` - Alignment-based extraction (primary method)
- `extract-rpf` - Pattern matching + HMM segmentation + optional STAR verification
- `decide-trim` - Auto-configure trimming by combining methods
- `check-cleanliness` - Categorized cleanliness checking
- `detect-adapter` - Specific adapter sequence detection
- `align-detect` - STAR alignment + feature detection

### 2. rRNA Removal

| Config ID | Method | Reference | Key Params |
|-----------|--------|-----------|------------|
| `rrna_none` | No filtering | - | - |
| `rrna_bowtie_ensembl` | Bowtie1 | Ensembl GTF-extracted rRNAs | `bowtie_rrna_mismatches`, `bowtie_rrna_k` |
| `rrna_bowtie_silva` | Bowtie1 | SILVA database | `bowtie_rrna_mismatches`, `bowtie_rrna_k` |
| `rrna_ribodetector` | RiboDetector | ML-based detection | `ribodetector_model` |

**RiboDetector Options:**
- CPU or GPU mode
- Different model sizes
- No reference required (ML-based)

### 3. Alignment Parameters

| Config ID | Mode | Mismatches | Soft-clip | Description |
|-----------|------|------------|-----------|-------------|
| `align_endtoend_mm0` | EndToEnd | 0 | No | Strictest - exact matches only |
| `align_endtoend_mm1` | EndToEnd | 1 | No | 1 mismatch allowed |
| `align_endtoend_mm2` | EndToEnd | 2 | No | 2 mismatches allowed |
| `align_endtoend_mm3` | EndToEnd | 3 | No | Default - 3 mismatches |
| `align_local_mm3` | Local | 3 | Yes | Allows soft-clipping |
| `align_local_mm3_sc5` | Local | 3 | 5' only | 5' soft-clip tolerance |
| `align_local_mm3_sc3` | Local | 3 | 3' only | 3' soft-clip tolerance |

**STAR Parameters to Test:**
- `--alignEndsType`: EndToEnd vs Local
- `--outFilterMismatchNmax`: 0, 1, 2, 3
- `--clip5pNbases`, `--clip3pNbases`: Explicit trimming
- `--alignSoftClipAtReferenceEnds`: Yes/No

### 4. Multi-mapper Handling

| Config ID | Max Loci | Description |
|-----------|----------|-------------|
| `multimap_1` | 1 | Unique mappers only |
| `multimap_5` | 5 | Up to 5 locations |
| `multimap_10` | 10 | Default - up to 10 locations |
| `multimap_20` | 20 | Permissive |
| `multimap_100` | 100 | Very permissive |

### 5. Read Length Filtering

| Config ID | Range | Description |
|-----------|-------|-------------|
| `length_none` | All | No filtering (default) |
| `length_strict` | 28-32nt | Typical monosome range |
| `length_relaxed` | 25-35nt | Broader range |
| `length_euk_short` | 26-30nt | Short for some organisms |
| `length_euk_long` | 29-33nt | Longer monosome |

### 6. P-site Offset Calculation

| Config ID | Tool | Method | Key Params |
|-----------|------|--------|------------|
| `offset_ribowaltz_auto` | RiboWaltz | Auto-select extremity | `flanking=6` |
| `offset_ribowaltz_5end` | RiboWaltz | Force 5' end | `extremity='5end'` |
| `offset_ribowaltz_3end` | RiboWaltz | Force 3' end | `extremity='3end'` |
| `offset_ribometric_changepoint` | RiboMetric | Changepoint detection | `--offset-calculation-method changepoint` |
| `offset_ribometric_ribowaltz` | RiboMetric | RiboWaltz-style | `--offset-calculation-method ribowaltz` |
| `offset_external_ribowaltz` | RiboMetric | Pass RiboWaltz offsets | `--offset-read-length` |
| `offset_global_15` | RiboMetric | Fixed offset | `--offset-global 15` |
| `offset_global_12` | RiboMetric | Fixed offset | `--offset-global 12` |

**RiboMetric Offset Options:**
- `--offset-calculation-method`: changepoint, ribowaltz, tripsviz
- `--offset-read-length`: TSV file with length-specific offsets
- `--offset-read-specific`: TSV file with read-specific offsets
- `--offset-global`: Single global offset value

---

## Metrics to Collect

### Stage 1: Pre-alignment Metrics

| Metric | Source | Purpose |
|--------|--------|---------|
| Raw read count | FASTQ | Starting point |
| Post-collapse count | Collapsed FASTA | Unique sequences |
| Adapter detection rate | getRPF/fastp | Adapter identification success |
| Detected adapter sequence | getRPF | What was found |
| Architecture confidence | getRPF | Extraction confidence |
| rRNA removal rate | Bowtie/RiboDetector | Contamination level |
| rRNA false positive rate | Test on CDS reads | Over-filtering? |

### Stage 2: Alignment Metrics

| Metric | Source | Purpose |
|--------|--------|---------|
| Total aligned reads | STAR Log.final.out | Mapping success |
| Uniquely mapped reads | STAR | Primary signal |
| Multi-mapped reads | STAR | Repetitive regions |
| Unmapped reads | STAR | Data loss |
| Mismatch distribution | BAM MD tag | Sequence quality |
| Soft-clip length distribution | CIGAR | Trimming adequacy |
| Alignment to CDS vs non-CDS | RiboMetric | Signal specificity |

### Stage 3: Quality Metrics (RiboMetric)

| Metric | Source | Good Value | Purpose |
|--------|--------|------------|---------|
| Frame 0 percentage | RiboMetric | >65% | Triplet periodicity |
| Frame 1 percentage | RiboMetric | <20% | |
| Frame 2 percentage | RiboMetric | <20% | |
| Information content | RiboMetric | >1.5 bits | Data richness |
| Read length peak | RiboMetric | 28-32nt | Expected RPF size |
| Start codon enrichment | RiboMetric | >3x background | Translation signal |
| Stop codon profile | RiboMetric | Sharp drop | Termination signal |
| CDS/UTR ratio | RiboMetric | >5:1 | Specificity |
| Metagene periodicity | RiboMetric | Clear 3nt pattern | Translation signature |

### Stage 4: Offset Quality (RiboWaltz)

| Metric | Source | Purpose |
|--------|--------|---------|
| Offset by length | RiboWaltz | Per-length precision |
| Offset confidence | RiboWaltz | Statistical robustness |
| Extremity selection | RiboWaltz | 5' vs 3' consistency |
| Correlation with start | RiboWaltz | Biological validation |

### Stage 5: Coverage Metrics

| Metric | Source | Purpose |
|--------|--------|---------|
| Per-gene coverage depth | BigWig | Signal strength |
| Coverage uniformity | BigWig | Bias detection |
| 5'UTR/CDS/3'UTR ratios | BEDgraph | Translation focus |
| Strand balance | BEDgraph | Library quality |

---

## Visualization Requirements

### 1. Read Flow Sankey Diagram
**Purpose:** Show data flow through pipeline stages
**Data needed:**
- Read counts at each stage
- Split by: kept/filtered/failed
- Color by fate

### 2. Metric Comparison Heatmaps
**Purpose:** Compare configs across multiple metrics
**Layout:**
- Rows: Configuration names
- Columns: Metrics (normalized)
- Color: Red (bad) → White → Green (good)

### 3. Frame Distribution Bar Charts
**Purpose:** Show triplet periodicity
**Layout:**
- Grouped bars by config
- Three bars per group (F0, F1, F2)
- Highlight >65% F0 threshold

### 4. Metagene Overlay Plots
**Purpose:** Compare start/stop codon profiles
**Layout:**
- Multiple lines (one per config)
- Shaded confidence intervals
- Vertical line at start/stop codon

### 5. Read Length Distributions
**Purpose:** Show RPF size selection
**Layout:**
- Overlaid histograms
- Highlight 28-32nt range
- Show peak for each config

### 6. Offset Comparison Scatter
**Purpose:** Compare offset methods
**Layout:**
- X: RiboWaltz offset
- Y: RiboMetric offset
- Points: read lengths
- Diagonal = agreement

### 7. Alignment Rate Comparison
**Purpose:** Show data retention across configs
**Layout:**
- Stacked bar chart
- Unique / Multi / Unmapped
- Per config

### 8. rRNA Contamination Plot
**Purpose:** Compare rRNA removal effectiveness
**Layout:**
- Before/after comparison
- By method
- False positive rate overlay

---

## Output File Structure

```
validation/
├── results/
│   ├── {config_name}/              # Per-config pipeline outputs
│   │   ├── collapsed/
│   │   ├── alignment/
│   │   ├── RiboMetric/
│   │   ├── ribowaltz/
│   │   └── bigwigs/
│   └── metrics/
│       ├── all_metrics.tsv         # Combined metrics table
│       ├── stage1_preAlignment.tsv
│       ├── stage2_alignment.tsv
│       ├── stage3_quality.tsv
│       ├── stage4_offset.tsv
│       └── stage5_coverage.tsv
├── figures/
│   ├── sankey_read_flow.svg
│   ├── heatmap_metrics.svg
│   ├── bars_frame_distribution.svg
│   ├── metagene_start.svg
│   ├── metagene_stop.svg
│   ├── hist_read_lengths.svg
│   ├── scatter_offsets.svg
│   ├── bars_alignment_rates.svg
│   └── rrna_comparison.svg
├── reports/
│   ├── validation_summary.html     # Interactive report
│   └── validation_slides.pdf       # Presentation-ready
└── configs/
    └── *.yaml                      # All test configs
```

---

## Recommended Test Datasets

### Minimal Test Set (fast iteration)
- 1 yeast sample (S. cerevisiae) - small genome, fast
- 1 human sample - real-world complexity

### Comprehensive Test Set
| Organism | Study | Samples | Library Type | Purpose |
|----------|-------|---------|--------------|---------|
| S. cerevisiae | SRP000001 | 3 | Standard | Baseline |
| H. sapiens | SRP000002 | 3 | Standard | Human reference |
| M. musculus | SRP000003 | 2 | Standard | Mammalian |
| H. sapiens | SRP000004 | 2 | Flash-freeze | Protocol variation |
| H. sapiens | SRP000005 | 2 | Harringtonine | Drug treatment |

*(Replace SRP numbers with actual accessions)*

---

## Decision Framework

### For each decision point, document:

1. **Options tested** (list configs)
2. **Key metrics compared** (table)
3. **Statistical test** (if applicable)
4. **Recommendation** (with justification)
5. **Caveats** (when recommendation might not apply)

### Example Decision: Alignment Mismatch Tolerance

```
DECISION: Use 2 mismatches for standard processing

EVIDENCE:
- 0 mismatches: Loses 15% of reads, no quality improvement
- 1 mismatch: Loses 8% of reads, minimal quality gain
- 2 mismatches: Optimal balance (98% retention, 68% F0)
- 3 mismatches: Marginal gain, slight decrease in F0

CAVEATS:
- For high-quality libraries, 1 mismatch may be preferred
- For degraded samples, 3 mismatches may recover more signal
```

---

## Implementation Priority

### Phase 1: Core Comparisons
1. Adapter detection methods (4 configs)
2. Alignment stringency (4 configs)
3. Multi-mapper handling (5 configs)

### Phase 2: Refinements
4. rRNA removal methods (4 configs)
5. Read length filtering (5 configs)
6. Offset calculation methods (7 configs)

### Phase 3: Advanced
7. Soft-clipping tolerance (3 configs)
8. Combined optimal configurations
9. Organism-specific optimization

---

## Next Steps

1. **Create all config YAML files** for each validation dimension
2. **Build metrics collection script** to aggregate outputs
3. **Create figure generation notebook** using plotly/matplotlib
4. **Run validation on test datasets**
5. **Generate report and slides**
