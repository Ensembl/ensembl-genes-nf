# Ensembl Statistics Workflow

The Ensembl Statistics workflow generates comprehensive gene set metrics for Ensembl core databases, producing standardized statistics used across all Ensembl projects.

## Overview

This workflow computes essential statistics about genes, transcripts, exons, and other genomic features stored in Ensembl core databases. These statistics are critical for:

- Database quality assessment
- Public data displays on Ensembl website
- Comparative genomics analyses
- Release validation and consistency checks

## Workflow Components

### 1. Statistics Generation

Computes core metrics from the database:

```bash
nextflow run main.nf \
  --csvFile databases.csv \
  --run_ensembl_stats \
  --host mysql-server.example.com \
  --port 3306 \
  --user_r ensro \
  --enscode /path/to/ensembl/modules \
  --outdir results
```

### 2. Beta Metakeys

Generates metadata for beta releases:

```bash
nextflow run main.nf \
  --csvFile databases.csv \
  --run_ensembl_beta_metakeys \
  --host mysql-server.example.com \
  --user_r ensro \
  --enscode /path/to/ensembl/modules
```

### 3. Apply to Database

Load generated statistics back into the database:

```bash
nextflow run main.nf \
  --csvFile databases.csv \
  --run_ensembl_stats \
  --apply_ensembl_stats \
  --host mysql-server.example.com \
  --user ensadmin \
  --password secret123 \
  --enscode /path/to/ensembl/modules \
  --team genebuild
```

## Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--enscode` | Path to Ensembl API modules | `/nfs/software/ensembl/ENSCODE` |
| `--host` | Database host server | `mysql-ens-sta-5.ebi.ac.uk` |
| `--port` | Database port | `4686` |
| `--user_r` | Read-only user (for generation) | `ensro` |
| `--user` | Write user (for applying) | `ensadmin` |
| `--password` | Database password (for applying) | - |
| `--team` | Team responsible (metakey) | `genebuild` |

### CSV Input

```csv
dbname,species_id
homo_sapiens_core_110_38,1
mus_musculus_core_110_39,1
danio_rerio_core_110_11,1
```

## Generated Statistics

### Gene-Level Metrics

| Statistic | Description | Example Value |
|-----------|-------------|---------------|
| `coding_cnt` | Number of protein-coding genes | 19,950 |
| `pseudogene_cnt` | Number of pseudogenes | 14,723 |
| `noncoding_cnt_*` | Non-coding genes by biotype | varies |
| `gene_cnt` | Total gene count | 60,670 |
| `alt_gene_cnt` | Genes on alternate sequences | 2,156 |

### Transcript-Level Metrics

| Statistic | Description | Example Value |
|-----------|-------------|---------------|
| `coding_transcript_cnt` | Protein-coding transcripts | 87,521 |
| `noncoding_transcript_cnt` | Non-coding transcripts | 51,210 |
| `transcript_cnt` | Total transcripts | 251,769 |
| `alt_transcript_cnt` | Transcripts on alternate seqs | 8,452 |

### Structural Metrics

| Statistic | Description | Example Value |
|-----------|-------------|---------------|
| `exon_cnt` | Total exons | 756,234 |
| `intron_cnt` | Total introns | 504,465 |
| `coding_exon_cnt` | Coding exons | 512,890 |

### Assembly Statistics

| Statistic | Description | Example Value |
|-----------|-------------|---------------|
| `toplevel_seq_cnt` | Top-level sequences | 639 |
| `chromosome_cnt` | Chromosomes | 24 |
| `scaffold_cnt` | Scaffolds | 615 |
| `contig_cnt` | Contigs | 50,543 |

### Gene Set Characteristics

| Statistic | Description | Example Value |
|-----------|-------------|---------------|
| `known_gene_cnt` | Genes with external refs | 18,234 |
| `novel_gene_cnt` | Novel predictions | 1,716 |
| `avg_transcript_per_gene` | Mean transcripts/gene | 4.15 |
| `avg_exon_per_transcript` | Mean exons/transcript | 8.7 |

## Output Files

### Statistics JSON

`{database}_statistics.json`

```json
{
  "database": "homo_sapiens_core_110_38",
  "species_id": 1,
  "generated_date": "2024-01-15T10:30:00",
  "statistics": {
    "coding_cnt": 19950,
    "pseudogene_cnt": 14723,
    "noncoding_cnt_lncRNA": 17910,
    "noncoding_cnt_misc_RNA": 2212,
    "noncoding_cnt_miRNA": 1879,
    "noncoding_cnt_rRNA": 549,
    "noncoding_cnt_snoRNA": 943,
    "noncoding_cnt_snRNA": 1901,
    "gene_cnt": 60670,
    "transcript_cnt": 251769,
    "coding_transcript_cnt": 87521,
    "exon_cnt": 756234,
    "toplevel_seq_cnt": 639,
    "chromosome_cnt": 24
  }
}
```

### Directory Structure

```
results/
└── ensembl_stats/
    ├── homo_sapiens_core_110_38_statistics.json
    ├── mus_musculus_core_110_39_statistics.json
    └── danio_rerio_core_110_11_statistics.json
```

## Use Cases

### 1. Release Validation

Validate statistics before public release:

```bash
# Generate statistics for all release databases
nextflow run main.nf \
  --csvFile release_110_databases.csv \
  --run_ensembl_stats \
  --host mysql-ens-sta-5.ebi.ac.uk \
  --port 4686 \
  --user_r ensro \
  --enscode /nfs/software/ensembl/ENSCODE \
  --outdir release_110_stats
```

**Check for:**
- Gene count consistency with previous release
- Reasonable transcript-to-gene ratios
- Expected biotype distributions
- Assembly sequence counts match expectations

### 2. Database Updates

Generate and apply statistics to updated databases:

```bash
nextflow run main.nf \
  --csvFile updated_databases.csv \
  --run_ensembl_stats \
  --apply_ensembl_stats \
  --host mysql-ens-sta-5.ebi.ac.uk \
  --port 4686 \
  --user ensadmin \
  --password ${DB_PASSWORD} \
  --enscode /nfs/software/ensembl/ENSCODE \
  --team genebuild \
  --project ensembl
```

### 3. Comparative Analysis

Compare statistics across species:

```bash
# Generate stats for a taxonomic group
cat > vertebrates.csv << EOF
dbname,species_id
homo_sapiens_core_110_38,1
mus_musculus_core_110_39,1
gallus_gallus_core_110_7,1
xenopus_tropicalis_core_110_10,1
danio_rerio_core_110_11,1
EOF

nextflow run main.nf \
  --csvFile vertebrates.csv \
  --run_ensembl_stats \
  --host mysql-server.example.com \
  --user_r ensro \
  --enscode /software/ensembl/ENSCODE
```

Then analyze:

```python
import json
import pandas as pd

# Load all statistics
stats = []
for f in glob.glob("results/ensembl_stats/*.json"):
    with open(f) as fh:
        stats.append(json.load(fh))

# Create comparison DataFrame
df = pd.DataFrame([s['statistics'] for s in stats])
df['species'] = [s['database'].split('_core_')[0] for s in stats]

# Compare metrics
print(df[['species', 'gene_cnt', 'coding_cnt', 'avg_transcript_per_gene']])
```

### 4. Quality Metrics Tracking

Track metrics over time:

```bash
# Generate statistics for version tracking
cat > versions.csv << EOF
dbname,species_id
species_v1_core,1
species_v2_core,1
species_v3_core,1
EOF

nextflow run main.nf \
  --csvFile versions.csv \
  --run_ensembl_stats \
  --host mysql-server.example.com \
  --user_r ensro \
  --enscode /software/ensembl/ENSCODE
```

## Beta Metakeys

### What are Beta Metakeys?

Metadata entries for beta/pre-release databases that include:

- Release preparation status
- Data freeze dates
- Validation checkpoints
- Team responsibilities

### Generating Beta Metakeys

```bash
nextflow run main.nf \
  --csvFile beta_databases.csv \
  --run_ensembl_beta_metakeys \
  --host mysql-ens-staging.ebi.ac.uk \
  --user_r ensro \
  --enscode /software/ensembl/ENSCODE \
  --team genebuild
```

### Applying Beta Metakeys

```bash
nextflow run main.nf \
  --csvFile beta_databases.csv \
  --run_ensembl_beta_metakeys \
  --apply_ensembl_beta_metakeys \
  --host mysql-ens-staging.ebi.ac.uk \
  --user ensadmin \
  --password ${DB_PASSWORD} \
  --enscode /software/ensembl/ENSCODE \
  --team genebuild \
  --project ensembl
```

## Interpreting Statistics

### Gene Count Expectations

| Organism Type | Typical Gene Count | Notes |
|---------------|-------------------|-------|
| Mammals | 20,000 - 25,000 | Fairly consistent |
| Birds | 15,000 - 20,000 | Slightly lower than mammals |
| Fish | 20,000 - 30,000 | Variable; some WGD events |
| Insects | 10,000 - 20,000 | *Drosophila* ~14k |
| Plants | 25,000 - 40,000 | Often polyploid |
| Fungi | 5,000 - 15,000 | Yeast ~6k |

### Transcript-to-Gene Ratios

| Ratio | Interpretation |
|-------|---------------|
| 1.0 - 2.0 | Simple organisms or minimal annotation |
| 2.0 - 5.0 | Typical for well-annotated genomes |
| 5.0 - 10.0 | Complex alternative splicing (e.g., human) |
| >10.0 | Possible over-prediction; review |

### Coding vs. Non-coding

| Organism | Typical Coding % | Notes |
|----------|-----------------|-------|
| Mammals | 30-40% | Many lncRNAs annotated |
| Fish | 40-60% | Variable annotation depth |
| Insects | 60-80% | Fewer non-coding annotated |
| Yeast | 80-90% | Compact genomes |

## Applying Statistics to Database

### Prerequisites

- Write access to database (`--user` and `--password`)
- Ensembl API modules configured (`--enscode`)
- Valid `meta` table in core database

### Process

1. **Generate statistics** from database queries
2. **Validate** computed metrics
3. **Insert** into `meta` table with appropriate keys
4. **Update** display names and descriptions

### Metakey Format

Statistics are stored as `meta` table entries:

| meta_key | meta_value |
|----------|------------|
| `genebuild.statistics.gene_cnt` | `60670` |
| `genebuild.statistics.coding_cnt` | `19950` |
| `genebuild.statistics.transcript_cnt` | `251769` |

### Team Attribution

The `--team` parameter adds responsibility metadata:

```sql
INSERT INTO meta (meta_key, meta_value) 
VALUES ('genebuild.team', 'genebuild');
```

## Troubleshooting

### Missing Ensembl API

**Problem:** Pipeline cannot find Ensembl modules

```
Error: Can't locate Bio/EnsEMBL/DBSQL/DBAdaptor.pm
```

**Solution:** Provide correct path to ENSCODE:

```bash
# Check Ensembl installation
ls /nfs/software/ensembl/ENSCODE/ensembl/modules

# Set in pipeline
--enscode /nfs/software/ensembl/ENSCODE
```

### Database Connection Failed

**Problem:** Cannot connect to database

**Solutions:**
```bash
# Test connection
mysql -h mysql-server.example.com -P 3306 -u ensro -e "SHOW DATABASES;"

# Check credentials
--host mysql-server.example.com \
--port 3306 \
--user_r ensro
```

### BioPerl Issues

**Problem:** BioPerl not found

**Solution:** Specify BioPerl path:

```bash
--bioperl /usr/local/bioperl-1.6.924
```

### Permission Denied (Apply Mode)

**Problem:** Cannot write to database

**Solution:** Use user with write permissions:

```bash
--user ensadmin \
--password ${DB_PASSWORD}
```

### Unexpected Statistics

**Problem:** Gene counts seem wrong

**Solutions:**
1. Check `species_id` is correct
2. Verify database schema version
3. Check for PAR regions (pseudoautosomal)
4. Review gene biotype classifications

## Best Practices

!!! tip "Generate Before Applying"
    Always run with `--run_ensembl_stats` alone first to review statistics before applying to the database.

!!! tip "Version Control Statistics"
    Keep JSON outputs in version control to track changes across releases.

!!! tip "Validate Against Previous Release"
    Compare statistics to previous release to catch unexpected changes:
    ```bash
    diff release_109_stats/ release_110_stats/
    ```

!!! tip "Use Read-Only User for Generation"
    Generate statistics with `--user_r` (read-only) to avoid accidental modifications.

!!! tip "Backup Before Applying"
    Always backup the database before applying statistics:
    ```bash
    mysqldump -h host -u user -p database > backup.sql
    ```

## Advanced Topics

### Custom Statistics Queries

Modify or add statistics by editing SQL query files:

```bash
--meta_query_file custom_queries.sql
```

### Batch Processing

Process hundreds of databases:

```bash
# Get all core databases for a release
mysql -h mysql-server.example.com -u ensro -e "
  SELECT schema_name, 1 
  FROM information_schema.schemata 
  WHERE schema_name LIKE '%_core_110_%'
" > release_110_cores.csv

# Add header
sed -i '1s/^/dbname,species_id\n/' release_110_cores.csv

# Run pipeline
nextflow run main.nf \
  --csvFile release_110_cores.csv \
  --run_ensembl_stats \
  --host mysql-server.example.com \
  --user_r ensro \
  --enscode /software/ensembl/ENSCODE
```

### Integration with Ensembl Production

For Ensembl production pipelines:

```bash
nextflow run main.nf \
  --csvFile ${ENSEMBL_RELEASE}/databases.csv \
  --run_ensembl_stats \
  --apply_ensembl_stats \
  --host ${ENSEMBL_STAGING_HOST} \
  --user ${ENSEMBL_ADMIN_USER} \
  --password ${ENSEMBL_ADMIN_PASS} \
  --enscode ${ENSCODE} \
  --team genebuild \
  --project ensembl \
  --mysql_ensadmin /software/ensembl/ensadmin
```

## Performance

### Resource Requirements

| Database Size | Memory | CPU | Runtime |
|---------------|--------|-----|---------|
| Small (<10k genes) | 4 GB | 2 | 1-2 min |
| Medium (10-30k genes) | 8 GB | 4 | 2-5 min |
| Large (>30k genes) | 16 GB | 8 | 5-10 min |

### Optimization

Statistics generation is database I/O bound:

- Use local database connections when possible
- Ensure database server has adequate resources
- Process multiple databases in parallel

## Next Steps

- [BUSCO Workflow](busco.md) - Gene set completeness assessment
- [OMArk Workflow](omark.md) - Proteome quality validation
- [Output Documentation](../output.md) - Complete output reference
- [Troubleshooting](../troubleshooting.md) - Common issues and solutions
