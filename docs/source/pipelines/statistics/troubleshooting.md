# Statistics Pipeline Troubleshooting Guide

This guide covers common issues, error messages, and solutions for the Ensembl genes statistics pipeline.

## Table of Contents

1. [Database Connection Issues](#database-connection-issues)
2. [Module-Specific Errors](#module-specific-errors)
3. [Resource and Performance Issues](#resource-and-performance-issues)
4. [File System and I/O Problems](#file-system-and-io-problems)
5. [Version and Dependency Issues](#version-and-dependency-issues)
6. [Data Quality and Validation](#data-quality-and-validation)
7. [Debugging Strategies](#debugging-strategies)

---

## Database Connection Issues

### Error: "Access denied for user"

**Symptoms:**
```
ERROR: Access denied for user 'username'@'host' (using password: YES)
```

**Causes:**
- Incorrect database credentials
- Missing database privileges
- Network access restrictions

**Solutions:**
1. Verify credentials in your configuration:
   ```groovy
   params.user = 'correct_username'
   params.password = 'correct_password'
   params.host = 'correct_host'
   params.port = 3306
   ```

2. Test connection manually:
   ```bash
   mysql -h ${params.host} -P ${params.port} -u ${params.user} -p
   ```

3. Check database grants:
   ```sql
   SHOW GRANTS FOR 'username'@'host';
   ```

**Affected Modules:** All database-dependent modules (FETCH_GENOME, FETCH_PROTEINS, RUN_STATISTICS, RUN_ENSEMBL_META, POPULATE_DB, BUSCO_CORE_METAKEYS)

---

### Error: "Can't connect to MySQL server"

**Symptoms:**
```
ERROR: Can't connect to MySQL server on 'host' (110)
```

**Causes:**
- Database server is down
- Network connectivity issues
- Firewall blocking connection
- Wrong host/port

**Solutions:**
1. Verify server is running:
   ```bash
   ping ${params.host}
   telnet ${params.host} ${params.port}
   ```

2. Check firewall rules allow outbound connections to database port

3. Verify host and port are correct in parameters

4. For compute environments, ensure security groups allow database access

---

### Error: "Unknown database"

**Symptoms:**
```
ERROR: Unknown database 'dbname'
```

**Causes:**
- Database name doesn't exist
- Typo in database name
- Metadata contains incorrect dbname

**Solutions:**
1. List available databases:
   ```sql
   SHOW DATABASES LIKE '%pattern%';
   ```

2. Verify metadata `dbname` field matches actual database name

3. Check for case sensitivity (database names are case-sensitive on some systems)

4. Ensure database was created before running pipeline

---

## Module-Specific Errors

### BUSCO Modules

#### Error: "BUSCO dataset not found"

**Symptoms:**
```
ERROR: Could not download dataset 'lineage_odb10'
```

**Causes:**
- Network connectivity to BUSCO database
- Invalid lineage name
- Corrupted cache directory

**Solutions:**
1. Verify lineage name is valid:
   ```bash
   busco --list-datasets
   ```

2. Clear cache and retry:
   ```bash
   rm -rf ${params.cacheDir}/busco_downloads/
   ```

3. Check internet connectivity:
   ```bash
   curl -I https://busco-data.ezlab.org/
   ```

4. Manually download dataset if needed:
   ```bash
   busco --download ${params.busco_lineage}
   ```

**Affected Modules:** BUSCO_DATASET, BUSCO_GENOME_LINEAGE, BUSCO_PROTEIN_LINEAGE

---

#### Error: "BUSCO failed to run"

**Symptoms:**
```
ERROR: BUSCO analysis failed for sample
```

**Causes:**
- Input file format issues
- Insufficient memory
- Corrupted input data

**Solutions:**
1. Validate input file format:
   ```bash
   # For genome
   head -n 20 genome.fa
   
   # For proteins
   head -n 20 proteins.fa
   ```

2. Check file is not empty:
   ```bash
   wc -l input.fa
   ```

3. Increase memory allocation in config:
   ```groovy
   process {
       withLabel: busco {
           memory = '16 GB'
       }
   }
   ```

4. Check BUSCO logs in work directory for detailed error

---

### OMAmer/OMark Modules

#### Error: "OMAmer database not found"

**Symptoms:**
```
ERROR: [Errno 2] No such file or directory: 'omamer_database'
```

**Causes:**
- `params.omamer_database` not set correctly
- Database file doesn't exist
- Insufficient permissions

**Solutions:**
1. Verify database path exists:
   ```bash
   ls -lh ${params.omamer_database}
   ```

2. Set correct path in parameters:
   ```groovy
   params.omamer_database = '/path/to/LUCA.h5'
   ```

3. Download database if missing:
   ```bash
   # See OMAmer documentation for download instructions
   ```

**Affected Modules:** OMAMER_HOG, OMARK

---

#### Error: "OMAmer search failed"

**Symptoms:**
```
ERROR: Search failed with exit code 1
```

**Causes:**
- Empty or invalid protein file
- Corrupted database
- Memory issues

**Solutions:**
1. Validate protein FASTA:
   ```bash
   grep -c "^>" proteins.fa  # Count sequences
   ```

2. Check database integrity:
   ```bash
   omamer validate-db ${params.omamer_database}
   ```

3. Increase maxForks if memory constrained:
   ```groovy
   process {
       withName: OMAMER_HOG {
           maxForks = 10  // Reduce from 15
       }
   }
   ```

---

### Fetch Modules

#### Error: "No translations found"

**Symptoms:**
```
WARNING: dump_translations.pl produced empty file
```

**Causes:**
- Database has no protein_coding genes
- Wrong database selected
- Database missing translation data

**Solutions:**
1. Check gene count:
   ```sql
   SELECT COUNT(*) FROM gene WHERE biotype = 'protein_coding';
   ```

2. Check translation table:
   ```sql
   SELECT COUNT(*) FROM translation;
   ```

3. Verify correct database in metadata:
   ```groovy
   meta.dbname = 'correct_core_db_name'
   ```

**Affected Modules:** FETCH_PROTEINS

---

#### Error: "Genome fetch failed"

**Symptoms:**
```
ERROR: Failed to fetch genome from database
```

**Causes:**
- Missing DNA sequences in database
- Insufficient memory for large genomes
- Database connection timeout

**Solutions:**
1. Verify DNA data exists:
   ```sql
   SELECT COUNT(*) FROM dna;
   ```

2. Increase timeout for large genomes:
   ```groovy
   process {
       withName: FETCH_GENOME {
           time = '6h'
           memory = '16 GB'
       }
   }
   ```

**Affected Modules:** FETCH_GENOME

---

### Statistics Modules

#### Error: "Statistics generation failed"

**Symptoms:**
```
ERROR: generate_species_homepage_stats.pl failed
```

**Causes:**
- Database schema issues
- Missing required tables
- Perl API version mismatch

**Solutions:**
1. Verify database schema version:
   ```sql
   SELECT meta_value FROM meta WHERE meta_key = 'schema_version';
   ```

2. Check required tables exist:
   ```sql
   SHOW TABLES LIKE '%gene%';
   SHOW TABLES LIKE '%transcript%';
   ```

3. Ensure Ensembl API matches schema:
   ```bash
   perl -MBio::EnsEMBL::Registry -e 'print $Bio::EnsEMBL::Registry::VERSION'
   ```

**Affected Modules:** RUN_STATISTICS, RUN_ENSEMBL_META

---

### Database Population

#### Error: "SQL execution failed"

**Symptoms:**
```
ERROR: mysql returned non-zero exit code
```

**Causes:**
- SQL syntax errors
- Missing table/column
- Insufficient privileges
- Duplicate key violations

**Solutions:**
1. Test SQL file manually:
   ```bash
   mysql -h ${params.host} -u ${params.user} -p ${meta.dbname} < file.sql
   ```

2. Check for errors in SQL:
   ```bash
   grep -i error file.sql
   ```

3. Verify INSERT/UPDATE privileges:
   ```sql
   SHOW GRANTS FOR CURRENT_USER();
   ```

4. Check for duplicate entries:
   ```sql
   SELECT meta_key, COUNT(*) FROM meta GROUP BY meta_key HAVING COUNT(*) > 1;
   ```

**Affected Modules:** POPULATE_DB, BUSCO_CORE_METAKEYS

---

## Resource and Performance Issues

### Error: "OutOfMemoryError" or process killed

**Symptoms:**
```
ERROR: Process killed (exit code 137)
OutOfMemoryError: Java heap space
```

**Causes:**
- Insufficient memory allocation
- Memory leak
- Processing very large files

**Solutions:**
1. Increase memory in configuration:
   ```groovy
   process {
       withLabel: busco {
           memory = { 8.GB * task.attempt }
       }
       withLabel: omamer {
           memory = { 16.GB * task.attempt }
       }
   }
   ```

2. Enable automatic retry with more memory:
   ```groovy
   process {
       errorStrategy = { task.exitStatus == 137 ? 'retry' : 'terminate' }
       maxRetries = 3
   }
   ```

3. Monitor memory usage:
   ```bash
   # Check work directory for .command.log
   tail -f work/xx/yyyy/.command.log
   ```

---

### Issue: Pipeline very slow / processes waiting

**Symptoms:**
- Many processes in "PENDING" state
- Low CPU usage
- Processes waiting hours to start

**Causes:**
- Too many parallel forks
- Database connection bottleneck
- I/O bottleneck

**Solutions:**
1. Reduce maxForks for database-heavy processes:
   ```groovy
   process {
       withLabel: fetch_file {
           maxForks = 10  // Reduce from 20
       }
       withName: OMAMER_HOG {
           maxForks = 5   // Reduce from 15
       }
   }
   ```

2. Stagger job submission:
   ```groovy
   executor {
       submitRateLimit = '10 sec'
   }
   ```

3. Check database connection pool limits

4. Monitor I/O with:
   ```bash
   iostat -x 5
   ```

---

### Issue: "Too many open files"

**Symptoms:**
```
ERROR: Too many open files
```

**Causes:**
- System file descriptor limit reached
- Too many parallel processes
- File handles not being released

**Solutions:**
1. Increase file descriptor limit:
   ```bash
   ulimit -n 4096
   ```

2. Set in your shell profile:
   ```bash
   echo "ulimit -n 4096" >> ~/.bashrc
   ```

3. Reduce maxForks globally:
   ```groovy
   executor.queueSize = 50
   ```

---

## File System and I/O Problems

### Error: "No space left on device"

**Symptoms:**
```
ERROR: No space left on device
```

**Causes:**
- Work directory full
- Cache directory full
- Output directory full

**Solutions:**
1. Check disk usage:
   ```bash
   df -h
   du -sh work/ outdir/ cacheDir/
   ```

2. Enable cleaning:
   ```groovy
   params.clean = true
   params.clean_work_dir = true
   ```

3. Use different storage for cache:
   ```groovy
   params.cacheDir = '/large/storage/path/'
   ```

4. Clean up old work directories:
   ```bash
   nextflow clean -f
   ```

---

### Issue: File system latency causing failures

**Symptoms:**
- Files not found immediately after creation
- Intermittent "No such file or directory" errors
- Process succeeds on retry

**Causes:**
- Distributed/NFS file system sync delay
- High I/O load

**Solutions:**
1. Increase latency delay:
   ```groovy
   params.files_latency = 5  // Increase from default 1
   ```

2. Use faster storage for work directory:
   ```groovy
   workDir = '/fast/local/storage/work'
   ```

3. Enable error retry:
   ```groovy
   process {
       errorStrategy = 'retry'
       maxRetries = 2
   }
   ```

---

### Error: "Permission denied"

**Symptoms:**
```
ERROR: Permission denied: '/path/to/file'
```

**Causes:**
- Insufficient file permissions
- Wrong user/group ownership
- Read-only file system

**Solutions:**
1. Check permissions:
   ```bash
   ls -la /path/to/file
   ```

2. Fix ownership:
   ```bash
   chown -R user:group /path/to/directory
   ```

3. Set correct permissions:
   ```bash
   chmod -R 755 /path/to/directory
   ```

4. Verify write access to output/cache directories

---

## Version and Dependency Issues

### Error: "Command not found"

**Symptoms:**
```
ERROR: busco: command not found
ERROR: omamer: command not found
```

**Causes:**
- Tool not installed in container
- Wrong container image
- PATH not set correctly

**Solutions:**
1. Verify container has the tool:
   ```bash
   docker run <container> which busco
   ```

2. Check container definition in config:
   ```groovy
   process {
       withLabel: busco {
           container = 'ezlabgva/busco:v5.4.7_cv1'
       }
   }
   ```

3. For conda environments:
   ```bash
   conda list | grep busco
   ```

---

### Error: Version incompatibility

**Symptoms:**
```
ERROR: This script requires BUSCO v5.x
ERROR: Ensembl API version mismatch
```

**Causes:**
- Outdated tool version
- Wrong Ensembl release
- Script-API version mismatch

**Solutions:**
1. Check tool versions in versions.yml outputs

2. Update container to correct version:
   ```groovy
   container = 'ezlabgva/busco:v5.4.7_cv1'  // Specific version
   ```

3. Match Ensembl API to database schema:
   ```bash
   # Database schema version
   mysql> SELECT meta_value FROM meta WHERE meta_key='schema_version';
   
   # Use matching API version
   git checkout release/110  # Match schema version
   ```

---

### Issue: Python/Perl module not found

**Symptoms:**
```
ERROR: Can't locate Bio/EnsEMBL/Registry.pm
ERROR: ModuleNotFoundError: No module named 'omamer'
```

**Causes:**
- Missing dependencies in environment
- Wrong Python/Perl environment active

**Solutions:**
1. For Perl modules:
   ```bash
   perl -MCPAN -e 'install Bio::EnsEMBL::Registry'
   ```

2. For Python modules:
   ```bash
   pip install omamer omark
   ```

3. Ensure correct environment in container:
   ```dockerfile
   RUN pip install omamer==0.3.3 omark==0.3.0
   ```

4. Check PERL5LIB and PYTHONPATH:
   ```bash
   echo $PERL5LIB
   echo $PYTHONPATH
   ```

---

## Data Quality and Validation

### Issue: Empty BUSCO results

**Symptoms:**
- BUSCO completes but reports 0% completeness
- No genes found

**Causes:**
- Wrong lineage selected
- Input data format issues
- Corrupted input file

**Solutions:**
1. Verify lineage is appropriate:
   ```bash
   busco --list-datasets
   # Choose lineage matching your species
   ```

2. Check input file format:
   ```bash
   # Should be valid FASTA
   grep "^>" input.fa | head
   ```

3. Try broader lineage:
   ```groovy
   params.busco_lineage = 'eukaryota_odb10'  // Start broad
   ```

4. Validate file integrity:
   ```bash
   md5sum input.fa
   ```

---

### Issue: OMark reports high contamination

**Symptoms:**
- OMark warns of possible contamination
- Many "unexpected" proteins

**Causes:**
- Actual contamination in assembly
- Wrong taxonomic placement
- Horizontal gene transfer

**Solutions:**
1. Review OMark detailed output in `omark_output/` directory

2. Check for known contaminants:
   ```bash
   grep -i "bacteria\|virus" omark_output/*.txt
   ```

3. Run contamination screening on assembly:
   ```bash
   # Use tools like blobtools, NCBI FCS-GX
   ```

4. If legitimate, document in metadata

---

### Issue: Statistics don't match expectations

**Symptoms:**
- Gene counts seem wrong
- Missing expected data in statistics

**Causes:**
- Database not fully populated
- Wrong database selected
- Filtering parameters

**Solutions:**
1. Manually check database:
   ```sql
   SELECT biotype, COUNT(*) FROM gene GROUP BY biotype;
   SELECT biotype, COUNT(*) FROM transcript GROUP BY biotype;
   ```

2. Verify using correct database:
   ```bash
   # Check meta.dbname matches intended database
   ```

3. Review statistics SQL files:
   ```bash
   cat core_statistics/*.sql
   ```

4. Compare with previous runs if available

---

## Debugging Strategies

### Strategy 1: Enable trace and debugging

```groovy
// nextflow.config
trace {
    enabled = true
    file = 'trace.txt'
}

dag {
    enabled = true
    file = 'dag.html'
}

report {
    enabled = true
    file = 'report.html'
}
```

### Strategy 2: Examine work directories

```bash
# Find failed process work directory
find work/ -name .exitcode -exec grep -l '1' {} \; | head -1 | xargs dirname

# View logs
cd <work_dir>
cat .command.sh    # Command that was run
cat .command.log   # stdout
cat .command.err   # stderr
cat .command.trace # Resource usage
```

### Strategy 3: Run process manually

```bash
# Navigate to work directory
cd <work_dir>

# Run command directly
bash .command.sh

# Or step through script line by line
```

### Strategy 4: Reduce scope for testing

```groovy
// Test with single sample
meta_ch = channel.of([
    [gca: 'GCA_000001405.15', dbname: 'test_db', production_name: 'test']
])
```

### Strategy 5: Check versions compatibility

```bash
# Collect all versions
cat work/**/versions.yml > all_versions.yml

# Compare against known working versions
```

### Strategy 6: Increase verbosity

```bash
# Run with debug logging
nextflow run main.nf -profile test --debug

# Nextflow trace
nextflow run main.nf -with-trace -with-report -with-timeline
```

---

## Common Error Patterns

### Pattern: Intermittent failures

**Indicators:**
- Process succeeds on retry
- Random timing of failures
- Different processes failing

**Likely Causes:**
- File system latency
- Network issues
- Resource contention

**Solution:**
```groovy
process {
    errorStrategy = 'retry'
    maxRetries = 2
}
params.files_latency = 5
```

---

### Pattern: All processes fail at same stage

**Indicators:**
- All samples fail at same module
- Consistent error message
- Happens immediately on start

**Likely Causes:**
- Configuration error
- Missing parameter
- Wrong path/credential

**Solution:**
- Review parameters for that specific module
- Check parameter documentation
- Validate paths and credentials

---

### Pattern: Failures after long runtime

**Indicators:**
- Process runs for hours then fails
- Memory-related errors
- Disk space errors

**Likely Causes:**
- Insufficient resources
- Memory leak
- Large file handling

**Solution:**
```groovy
process {
    memory = { 8.GB * task.attempt }
    time = { 4.h * task.attempt }
    errorStrategy = 'retry'
    maxRetries = 3
}
```

---

## Getting Help

### Before asking for help:

1. ✅ Check this troubleshooting guide
2. ✅ Review module documentation
3. ✅ Examine work directory logs
4. ✅ Check parameter configuration
5. ✅ Try with single sample/minimal data

### When reporting issues, include:

1. **Error message** (full text from .command.err)
2. **Command executed** (.command.sh contents)
3. **Configuration** (relevant params)
4. **Environment** (Nextflow version, executor)
5. **Steps to reproduce**
6. **Versions file** (if generated)

### Useful diagnostic commands:

```bash
# Nextflow version
nextflow -version

# Process logs
find work/ -name .command.err -exec grep -l ERROR {} \;

# Resource usage
cat work/**/.command.trace

# Configuration used
nextflow config -profile <your_profile>

# List failed processes
nextflow log <run_name> -f status,name,exit | grep FAILED
```

---

## Quick Reference

| Issue | First Check | Module |
|-------|-------------|---------|
| Can't connect to database | Credentials, network | All DB modules |
| BUSCO fails | Lineage, input format | BUSCO_* |
| OMAmer database not found | Path, file exists | OMAMER_HOG, OMARK |
| Empty output | Input data, database content | FETCH_*, BUSCO_* |
| SQL errors | Privileges, syntax | POPULATE_DB |
| Out of memory | Process memory config | BUSCO_*, OMAMER_* |
| File not found | File system latency | All |
| Slow pipeline | maxForks, database load | All |

---

**Last Updated**: 2026-02-07  
**For**: Ensembl Genes Statistics Pipeline v1.0