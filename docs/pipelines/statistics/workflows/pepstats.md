# Pepstats Workflow

The Pepstats workflow runs pepstats on ensembl protein sequences.

## Overview

The workflow computes essential statistics about protein sequences from Ensembl core databases.

## Workflow Components

Computes core metrics from the database:

```bash
nextflow run main.nf \
  --csvFile databases.csv \
  --run_pepstats \
  --host mysql-server.example.com \
  --port 3306 \
  --user_r ensro \
  --enscode /path/to/ensembl/modules \
  --outdir results
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

### Pepstats applied to the database

| Statistic | Description | Example Value |
|-----------|-------------|---------------|
| `IsoPoint` | Isoelectric point | 10.0024 |
| `Charge` | Charge | -0.5 |
| `MolecularWeight` | Number of residues | 100008.72 |
| `AvgResWeight` | Ave. residue weight | 100.092 |
| `NumResidues` | Number of residues | 100 |



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

### Permission Denied

**Problem:** Cannot write to database

**Solution:** Use user with write permissions:

```bash
--user ensadmin \
--password ${DB_PASSWORD}
```
