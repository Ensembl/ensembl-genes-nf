# BUSCO Workflow

**BUSCO** (Benchmarking Universal Single-Copy Orthologs) provides quantitative assessment of genome assembly and annotation completeness based on evolutionarily-informed expectations.

## Overview

BUSCO assesses genome quality by searching for conserved orthologous genes that are expected to be present in single copy in a given lineage. The presence, absence, or duplication of these genes provides insights into assembly and annotation quality.

## How BUSCO Works

1. **Lineage Selection**: Choose an appropriate lineage dataset based on taxonomy
2. **Gene Search**: Search for conserved single-copy orthologs using HMMER
3. **Classification**: Genes are classified as Complete, Duplicated, Fragmented, or Missing
4. **Scoring**: Generate completeness scores and detailed reports

## Analysis Modes

### Protein Mode

Analyzes translated protein sequences from gene predictions.

**When to use:**
- Assessing annotation quality
- Evaluating gene prediction completeness
- Comparing gene sets across species

**Input:** Protein FASTA from canonical transcripts

**Advantages:**
- Direct assessment of annotation quality
- Faster than genome mode
- Better for comparing gene sets

### Genome Mode

Analyzes the genome assembly directly using AUGUSTUS for gene prediction.

**When to use:**
- Assessing assembly completeness
- Checking for missing genes in annotation
- Validating genome assembly quality

**Input:** Genome FASTA (DNA sequence)

**Advantages:**
- Independent of existing annotation
- Can detect genes missed by annotation pipeline
- Better for assembly quality assessment

### Both Modes

Run both protein and genome modes for comprehensive assessment.

**When to use:**
- Complete quality control
- Comparing annotation vs. assembly completeness
- Production annotation validation

**Provides:**
- Annotation completeness (protein mode)
- Assembly completeness (genome mode)
- Comparison of both metrics

## Running BUSCO

### Using Core Database

```bash
# Protein mode only
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --busco_mode protein \
  --host mysql-server.example.com \
  --user_r ensro

# Genome mode only
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --busco_mode genome \
  --host mysql-server.example.com \
  --user_r ensro

# Both modes (recommended for QC)
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --busco_mode both \
  --host mysql-server.example.com \
  --user_r ensro
```

### Using NCBI Assembly

```bash
# Automatically download and analyze genome from NCBI
nextflow run main.nf \
  --csvFile ncbi_genomes.csv \
  --run_busco_ncbi \
  --outdir results
```

!!! info "NCBI Mode"
    NCBI mode only supports genome mode analysis. The assembly is automatically downloaded using the NCBI Datasets API.

## Choosing a Lineage

Select the **most specific** lineage that applies to your organism. More specific lineages provide better sensitivity.

### Decision Tree

```
Is your organism a...
├─ Eukaryote?
│  ├─ Animal (Metazoa)?
│  │  ├─ Vertebrate?
│  │  │  ├─ Mammal?
│  │  │  │  ├─ Primate? → primates_odb12
│  │  │  │  └─ Other → mammalia_odb12
│  │  │  ├─ Bird? → aves_odb12
│  │  │  ├─ Fish? → actinopterygii_odb12
│  │  │  └─ Other → vertebrata_odb12
│  │  └─ Invertebrate?
│  │     ├─ Insect?
│  │     │  ├─ Fly/Mosquito? → diptera_odb12
│  │     │  └─ Other → insecta_odb12
│  │     └─ Other → metazoa_odb12
│  ├─ Plant (Viridiplantae)?
│  │  ├─ Land plant? → embryophyta_odb12
│  │  └─ Other → viridiplantae_odb12
│  └─ Fungus?
│     ├─ Yeast? → saccharomycetes_odb12
│     └─ Other → fungi_odb12
└─ Bacteria? → bacteria_odb10
```

### Common Lineages

| Lineage | # BUSCOs | Description | Example Species |
|---------|----------|-------------|-----------------|
| `primates_odb12` | 13,780 | Primates | Human, chimp, gorilla |
| `mammalia_odb12` | 9,226 | Mammals | Mouse, dog, cow, whale |
| `aves_odb12` | 8,338 | Birds | Chicken, zebra finch |
| `actinopterygii_odb12` | 3,640 | Ray-finned fish | Zebrafish, medaka |
| `vertebrata_odb12` | 3,354 | Vertebrates | Any vertebrate |
| `diptera_odb12` | 3,285 | Flies and mosquitoes | *Drosophila*, *Anopheles* |
| `insecta_odb12` | 1,367 | Insects | Bee, beetle, butterfly |
| `embryophyta_odb12` | 1,614 | Land plants | Arabidopsis, rice |
| `fungi_odb12` | 758 | Fungi | Yeast, *Aspergillus* |
| `metazoa_odb12` | 954 | Animals | Any animal |
| `eukaryota_odb12` | 255 | Eukaryotes | Any eukaryote |

## Understanding Results

### BUSCO Categories

BUSCO classifies each searched gene into one of four categories:

| Category | Code | Description | Interpretation |
|----------|------|-------------|----------------|
| **Complete (C)** | C | Gene found in single copy with expected length | ✅ Good |
| **Complete and Duplicated (D)** | D | Gene found multiple times | ⚠️ May indicate duplication or misassembly |
| **Fragmented (F)** | F | Gene found but incomplete | ⚠️ Partial gene or assembly issue |
| **Missing (M)** | M | Gene not found | ❌ Gap in assembly/annotation |

### Score Interpretation

**Completeness Score = (C + D) / Total BUSCOs × 100**

| Score Range | Quality | Interpretation |
|-------------|---------|----------------|
| **95-100%** | Excellent | High-quality assembly/annotation |
| **90-95%** | Good | Acceptable for most analyses |
| **80-90%** | Fair | Some missing genes; investigate |
| **<80%** | Poor | Significant issues; troubleshoot |

!!! warning "Context Matters"
    Scores should be interpreted in context:
    
    - **Protein mode**: Lower scores may indicate annotation issues
    - **Genome mode**: Lower scores suggest assembly gaps
    - **Duplications**: High duplication (>5%) may indicate misassembly or polyploidy

### Example Output

#### Short Summary (`*_busco_short_summary.txt`)

```
# BUSCO version is: 5.4.7 
# The lineage dataset is: primates_odb12 (Creation date: 2024-01-08, number of genomes: 44, number of BUSCOs: 13780)
# Summarized benchmarking in BUSCO notation for file /data/proteins.fa
# BUSCO was run in mode: protein
# Gene predictor used: None

***** Results: *****

C:98.5%[S:97.8%,D:0.7%],F:0.8%,M:0.7%,n:13780	   
13570	Complete BUSCOs (C)			   
13474	Complete and single-copy BUSCOs (S)	   
96	Complete and duplicated BUSCOs (D)	   
110	Fragmented BUSCOs (F)			   
100	Missing BUSCOs (M)			   
13780	Total BUSCO groups searched		   
```

**Interpretation:**
- ✅ Excellent completeness (98.5%)
- ✅ Low duplication (0.7%)
- ✅ Very few fragmented (0.8%) or missing (0.7%)
- **Conclusion**: High-quality annotation

#### Full Table (`*_busco_full_table.tsv`)

```
# Busco id	Status	Sequence	Gene Start	Gene End	Strand	Score	Length
1000092at7742	Complete	ENSP00000262735	-	-	-	1256.4	416
1000093at7742	Complete	ENSP00000344818	-	-	-	2145.7	712
1000094at7742	Missing	-	-	-	-	-	-
1000095at7742	Fragmented	ENSP00000123456	-	-	-	456.2	189
```

## Output Files

### Per-Sample Outputs

| File | Description | Size |
|------|-------------|------|
| `{sample}_busco_short_summary.txt` | Summary statistics (protein mode) | ~1 KB |
| `{sample}_genome_busco_short_summary.txt` | Summary statistics (genome mode) | ~1 KB |
| `{sample}_busco_full_table.tsv` | Detailed results per BUSCO | ~1-5 MB |
| `{sample}_busco_missing_busco_list.tsv` | List of missing BUSCOs | Variable |

### Directory Structure

```
results/
└── busco/
    ├── sample1_busco_short_summary.txt
    ├── sample1_genome_busco_short_summary.txt
    ├── sample1_busco_full_table.tsv
    ├── sample2_busco_short_summary.txt
    └── sample2_genome_busco_short_summary.txt
```

## Interpreting Results by Mode

### Protein vs. Genome Comparison

| Scenario | Protein Score | Genome Score | Interpretation |
|----------|---------------|--------------|----------------|
| **Ideal** | High (>95%) | High (>95%) | Excellent assembly and annotation |
| **Good annotation** | High (>95%) | Medium (85-95%) | Good annotation, some assembly gaps |
| **Annotation issues** | Medium (85-95%) | High (>95%) | Assembly good, annotation needs improvement |
| **Both medium** | 85-95% | 85-95% | Both assembly and annotation have room for improvement |
| **Assembly problems** | Low (<85%) | Low (<85%) | Significant assembly gaps affecting annotation |

### High Duplication (>5%)

Possible causes:
- **Polyploidy**: Natural in some species (e.g., wheat, salmon)
- **Recent whole-genome duplication**: Normal in some lineages
- **Misassembly**: Collapsed repeats or chimeric contigs
- **Contamination**: Mixed samples
- **Incorrect lineage**: Using too broad a lineage

**Action:** Check genome characteristics and assembly quality metrics.

## Troubleshooting

### Low Completeness Scores

**Problem:** BUSCO completeness <80%

**Possible causes:**
1. **Wrong lineage**: Using inappropriate lineage dataset
2. **Assembly gaps**: Incomplete genome assembly
3. **Annotation issues**: Missing or incomplete gene predictions
4. **Contamination**: Non-target sequences in assembly
5. **Unusual biology**: Genuine gene loss (rare)

**Solutions:**
- Try a broader lineage (e.g., vertebrata instead of primates)
- Check assembly statistics (N50, gaps, contiguity)
- Review annotation pipeline parameters
- Screen for contamination

### High Fragmentation

**Problem:** Many fragmented BUSCOs (>5%)

**Possible causes:**
- Assembly fragmentation (low N50)
- Incorrect gene models
- Pseudogenes or processed pseudogenes

**Solutions:**
- Improve assembly contiguity
- Adjust gene prediction parameters
- Check for frame shifts or premature stop codons

### High Duplication

**Problem:** >10% duplicated BUSCOs (when not expected)

**Possible causes:**
- Haplotigs not removed during assembly
- Collapsed repeats
- Chimeric sequences
- Recent gene duplications (normal)

**Solutions:**
- Run purge_dups or similar tools
- Check assembly for haplotigs
- Validate with Hi-C or genetic maps

## Best Practices

!!! tip "Choose the right lineage"
    Use the **most specific lineage** available for your organism. More specific = better sensitivity and specificity.

!!! tip "Run both modes"
    For production annotations, always run both protein and genome modes to get complete picture of quality.

!!! tip "Compare across releases"
    Track BUSCO scores across annotation updates to monitor improvements or catch regressions.

!!! tip "Use with other QC metrics"
    BUSCO is one metric—combine with N50, BUSCO, gene count trends, and manual inspection.

!!! tip "Document exceptions"
    Some genuine biological variations (gene loss, lineage-specific duplications) can affect scores—document these.

## Advanced Topics

### Custom Lineages

For specialized projects, you can create custom lineage datasets:

1. Download from [OrthoDB](https://www.orthodb.org/)
2. Place in the lineage directory
3. Reference in your CSV or with `--busco_dataset`

### Offline Mode

Pre-download lineage datasets:

```bash
# Download lineages
busco --download_path /data/busco_lineages --download vertebrata_odb12

# Use in pipeline
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --download_path /data/busco_lineages
```

### Metakey Integration

Load BUSCO results as database metakeys:

```bash
nextflow run main.nf \
  --csvFile genomes.csv \
  --run_busco_core \
  --apply_busco_metakeys \
  --host mysql-server.example.com \
  --user ensadmin \
  --password secret123
```

## References

### Primary Citations

- **Manni et al. (2021).** BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. *Molecular Biology and Evolution*, 38(10):4647–4654. [DOI: 10.1093/molbev/msab199](https://doi.org/10.1093/molbev/msab199)

- **Manni et al. (2021).** BUSCO: Assessing genomic data quality and beyond. *Current Protocols*, 1:e323. [DOI: 10.1002/cpz1.323](https://doi.org/10.1002/cpz1.323)

### Additional Resources

- [BUSCO Website](https://busco.ezlab.org/)
- [BUSCO User Guide](https://busco.ezlab.org/busco_userguide.html)
- [OrthoDB Database](https://www.orthodb.org/)
- [BUSCO Publications](https://busco.ezlab.org/publications.html)

## Next Steps

- [OMArk Workflow](omark.md) - Complementary proteome assessment
- [Output Documentation](../output.md) - Detailed output file specifications
- [Troubleshooting](../troubleshooting.md) - Solve common issues
