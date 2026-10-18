# OMArk Workflow

**OMArk** (Orthology MARKer) is a proteome quality assessment tool that evaluates completeness based on the presence of conserved genes from the target lineage.

## Overview

OMArk provides an orthology-based assessment of proteome quality by comparing your protein set against a database of conserved orthologous groups. It identifies complete, fragmented, missing, and potentially contaminant sequences.

## How OMArk Works

1. **OMAmer Mapping**: Maps protein sequences to hierarchical orthologous groups (HOGs)
2. **Lineage Assignment**: Determines the taxonomic placement of each protein
3. **Completeness Assessment**: Evaluates expected gene presence based on lineage
4. **Consistency Check**: Identifies potential contamination or misassignments

## Key Features

- **Lineage-Specific**: Uses taxonomically appropriate conserved gene sets
- **Contamination Detection**: Identifies sequences inconsistent with target lineage
- **Fragmentation Analysis**: Distinguishes complete vs. partial genes
- **Placement Verification**: Validates taxonomic assignments

## Running OMArk

### Basic Usage

```bash
nextflow run main.nf \
  --csvFile proteomes.csv \
  --run_omark \
  --host mysql-server.example.com \
  --port 3306 \
  --user_r ensro \
  --outdir results
```

### Required Inputs

Create a CSV file with these columns:

```csv
dbname,species_id,taxon_id
danio_rerio_core_110_11,1,7955
xenopus_tropicalis_core_110_10,1,8364
gallus_gallus_core_110_7,1,9031
```

### Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--omamer_database` | System path | Path to OMAmer HOG database (.h5 file) |
| `--omark_singularity_path` | System path | Path to OMArk Singularity container |

### Custom Database

```bash
nextflow run main.nf \
  --csvFile proteomes.csv \
  --run_omark \
  --omamer_database /data/omamer/LUCA.h5 \
  --host mysql-server.example.com \
  --user_r ensro
```

## Understanding Results

### Output File

OMArk generates a detailed summary file: `{sample}_omark_proteins_detailed_summary.txt`

### Example Output

```
# OMArk Analysis Summary
# Database: danio_rerio_core_110_11
# Taxon ID: 7955 (Danio rerio)
# Date: 2024-01-15

Completeness Assessment:
========================
Expected HOGs: 3,640
Complete: 3,521 (96.7%)
Fragmented: 87 (2.4%)
Missing: 32 (0.9%)

Consistency Analysis:
====================
Consistent: 3,598 (98.8%)
Inconsistent: 10 (0.3%)
```

### Result Categories

| Category | Description | Interpretation |
|----------|-------------|----------------|
| **Complete** | Full-length conserved genes present | ✅ Expected genes found |
| **Fragmented** | Partial or truncated gene matches | ⚠️ Incomplete genes or pseudogenes |
| **Missing** | Expected genes not detected | ❌ Annotation gaps or genuine loss |
| **Consistent** | Proteins match expected lineage | ✅ Correct taxonomic placement |
| **Inconsistent** | Proteins from unexpected lineage | ⚠️ Contamination or HGT |

### Completeness Scores

**Overall Completeness = Complete / Expected × 100**

| Score Range | Quality | Interpretation |
|-------------|---------|----------------|
| **95-100%** | Excellent | High-quality, complete proteome |
| **90-95%** | Good | Minor gaps, generally good |
| **85-90%** | Fair | Some missing genes; investigate |
| **<85%** | Poor | Significant issues; troubleshoot |

### Consistency Scores

**Consistency = Consistent / Total × 100**

| Score Range | Quality | Interpretation |
|-------------|---------|----------------|
| **98-100%** | Excellent | Clean, no contamination |
| **95-98%** | Good | Minimal contamination |
| **90-95%** | Fair | Some inconsistent sequences |
| **<90%** | Poor | Significant contamination or misassignment |

## OMArk vs. BUSCO

Both tools assess proteome completeness but use different approaches:

| Feature | OMArk | BUSCO |
|---------|-------|-------|
| **Database** | OMAmer HOGs | OrthoDB single-copy orthologs |
| **Approach** | Hierarchical orthology | Single-copy gene presence |
| **Contamination Detection** | ✅ Built-in | ❌ Limited |
| **Lineage Specificity** | High | High |
| **Speed** | Fast | Slower |
| **Gene Count** | Variable by lineage | Fixed per lineage |
| **Best For** | Contamination screening | Assembly/annotation QC |

!!! tip "Use Both"
    OMArk and BUSCO complement each other. Use both for comprehensive quality assessment.

## Interpreting Results

### High Completeness, High Consistency

```
Complete: 97.5%
Consistent: 99.2%
```

**✅ Excellent Quality**
- Proteome is complete and clean
- No significant contamination
- Ready for downstream analyses

### High Completeness, Low Consistency

```
Complete: 96.8%
Consistent: 89.5%
```

**⚠️ Contamination Issues**
- Many genes present but from wrong lineage
- Possible contamination or mixed sample
- Check for foreign sequences

**Action:**
- Review inconsistent proteins
- Screen for contamination
- Verify sample purity

### Low Completeness, High Consistency

```
Complete: 82.3%
Consistent: 98.9%
```

**⚠️ Annotation Gaps**
- Correct lineage but missing genes
- Annotation incomplete or gene loss
- Assembly may have gaps

**Action:**
- Review annotation pipeline
- Check assembly completeness
- Consider re-annotation

### Low Completeness, Low Consistency

```
Complete: 78.5%
Consistent: 87.2%
```

**❌ Multiple Issues**
- Both completeness and contamination problems
- Severe quality issues
- Requires thorough investigation

**Action:**
- Comprehensive QC review
- Check sample identity
- Validate assembly and annotation

## Use Cases

### 1. Contamination Screening

OMArk excels at detecting contamination:

```bash
# Screen multiple assemblies
nextflow run main.nf \
  --csvFile new_assemblies.csv \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro
```

**Look for:**
- Low consistency scores (<95%)
- Inconsistent proteins from distantly related lineages
- Patterns suggesting systematic contamination

### 2. Pre-Publication QC

Validate proteome quality before public release:

```bash
# Run both OMArk and BUSCO
nextflow run main.nf \
  --csvFile release_genomes.csv \
  --run_omark \
  --run_busco_core \
  --busco_mode protein \
  --host mysql-server.example.com \
  --user_r ensro
```

### 3. Comparative Proteomics

Assess relative quality across multiple species:

```csv
dbname,species_id,taxon_id
species_a_core_110_1,1,12345
species_b_core_110_1,1,12346
species_c_core_110_1,1,12347
```

```bash
nextflow run main.nf \
  --csvFile comparative_set.csv \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro
```

### 4. Annotation Pipeline Validation

Test annotation pipeline improvements:

```bash
# Compare v1 vs. v2 annotations
cat > comparison.csv << EOF
dbname,species_id,taxon_id
species_v1_core,1,9606
species_v2_core,1,9606
EOF

nextflow run main.nf \
  --csvFile comparison.csv \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro
```

## Output Files

### Directory Structure

```
results/
└── omark/
    ├── sample1_omark_proteins_detailed_summary.txt
    ├── sample2_omark_proteins_detailed_summary.txt
    └── sample3_omark_proteins_detailed_summary.txt
```

### File Contents

The detailed summary includes:

1. **Header**: Metadata (database, taxon, date)
2. **Completeness Metrics**: HOG presence/absence
3. **Consistency Metrics**: Lineage validation
4. **Fragmentation Details**: Partial gene information
5. **Inconsistent Proteins**: List of potentially contaminating sequences

## Troubleshooting

### Low Completeness with BUSCO OK

**Problem:** OMArk shows low completeness but BUSCO is good

**Possible causes:**
- Different gene sets used (OMArk uses broader HOGs)
- Lineage-specific gene expansions/losses
- Different sensitivity thresholds

**Solution:** Review both results; some divergence is normal.

### High Inconsistency Rate

**Problem:** >10% inconsistent proteins

**Possible causes:**
1. **Contamination**: Foreign DNA in sample
2. **Wrong taxon ID**: Incorrect species assignment
3. **Horizontal gene transfer**: Natural (e.g., in bacteria)
4. **Symbionts**: Co-sequenced organisms

**Solutions:**
- Check sample purity
- Verify taxon ID is correct
- Screen for contamination with tools like BlobTools
- Review inconsistent sequences manually

### Missing OMAmer Database

**Problem:** Pipeline cannot find OMAmer database

**Solution:** Download and specify:

```bash
# Download from OMArk project
wget https://omabrowser.org/All/LUCA.h5

# Use in pipeline
nextflow run main.nf \
  --csvFile proteomes.csv \
  --run_omark \
  --omamer_database /data/LUCA.h5
```

### Very Long Runtime

**Problem:** OMArk taking too long

**Possible causes:**
- Very large proteome (>50,000 proteins)
- Slow I/O on database file

**Solutions:**
- Use local SSD for OMAmer database
- Split large proteomes if possible
- Increase computational resources

## Best Practices

!!! tip "Check Consistency First"
    Always review consistency scores before interpreting completeness. Contamination can inflate completeness scores.

!!! tip "Combine with BUSCO"
    Use OMArk and BUSCO together for comprehensive QC. They catch different types of issues.

!!! tip "Track Across Releases"
    Monitor OMArk scores across annotation updates to catch regressions or improvements.

!!! tip "Investigate Inconsistencies"
    Don't ignore inconsistent proteins—they often reveal real contamination or sample issues.

!!! tip "Use Appropriate Taxon IDs"
    Ensure taxon_id matches your species. Wrong IDs lead to incorrect consistency assessments.

## Advanced Topics

### Custom OMAmer Databases

For specialized analyses, you can create custom HOG databases:

1. Extract relevant clade from OMA
2. Build OMAmer index
3. Point pipeline to custom database

### Integration with Contamination Tools

Combine OMArk with other contamination detection:

```bash
# 1. Run OMArk
nextflow run main.nf --csvFile genomes.csv --run_omark

# 2. Extract inconsistent proteins
grep "Inconsistent" results/omark/*_detailed_summary.txt

# 3. Further screen with BlobTools or similar
```

### Batch Analysis

Process hundreds of proteomes:

```bash
# Generate CSV for all databases
mysql -h mysql-server.example.com -u ensro -e "
  SELECT CONCAT(db_name, ',1,', meta_value) 
  FROM information_schema.schemata 
  JOIN meta ON db_name = database()
  WHERE db_name LIKE '%_core_%' 
  AND meta_key = 'species.taxonomy_id'
" > all_cores.csv

# Run OMArk on all
nextflow run main.nf \
  --csvFile all_cores.csv \
  --run_omark \
  --host mysql-server.example.com \
  --user_r ensro
```

## Performance Considerations

### Resource Requirements

| Proteome Size | Memory | CPU | Runtime |
|---------------|--------|-----|---------|
| Small (<20k) | 8 GB | 4 | 10-20 min |
| Medium (20-40k) | 16 GB | 8 | 20-40 min |
| Large (>40k) | 32 GB | 16 | 40-90 min |

### Optimization

```bash
# Increase parallelization
nextflow run main.nf \
  --csvFile proteomes.csv \
  --run_omark \
  --max_cpus 32 \
  --max_memory 64.GB
```

## References

### Primary Citation

- **Nevers et al. (2022).** OMArk: proteome quality assessment using the OMA database. *bioRxiv*. [DOI: 10.1101/2022.11.25.517970](https://doi.org/10.1101/2022.11.25.517970)

### Additional Resources

- [OMArk GitHub](https://github.com/DessimozLab/OMArk)
- [OMAmer Documentation](https://omamer.readthedocs.io/)
- [OMA Browser](https://omabrowser.org/)

## Next Steps

- [BUSCO Workflow](busco.md) - Complementary completeness assessment
- [Ensembl Stats](ensembl-stats.md) - Generate database statistics
- [Output Documentation](../output.md) - Detailed output specifications
- [Troubleshooting](../troubleshooting.md) - Common issues and solutions
