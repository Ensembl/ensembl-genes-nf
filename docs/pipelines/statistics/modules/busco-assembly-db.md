# BUSCO Assembly DB

## Overview

The `BUSCO_ASSEMBLY_DB` process loads BUSCO **genome mode** results into the assembly metadata database. It parses the BUSCO short summary file and inserts or updates rows in the `assembly_metrics` table for the matching assembly, then sets `genome_busco.status` to `done` in the `assembly_events` table.

## Process Details

- **Label**: `python`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Cache**: Disabled because the process writes to the database
- **Publish Directory**: `${params.outdir}/${meta.gca}`
- **Conditional Execution**: Only runs when `params.update_registry` is true, and only on genome mode BUSCO results
- **Not with suppressed accessions**: The pipeline stops if `params.allow_suppressed_accessions` is also set. FETCH_GENOME may then download a suppressed GCA's `currentAccession` while `meta.gca` keeps the requested one, so the replacement's BUSCO would be recorded against the requested GCA

## Inputs

| Name | Type | Description |
|------|------|-------------|
| meta | val | Metadata map containing genome assembly info |
| summary_file | path | BUSCO genome summary file to parse |
| asm_metadata | val | JSON string with the assembly metadata DB connection params |

### Required Metadata Fields

- `gca`: Versioned genome assembly accession (e.g. `GCA_000001405.29`)
- `busco_mode`: Must be `genome`

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| metrics_json | tuple val(meta), path | `<gca>_busco_genome_asmdb.json` with the metrics loaded, published to `${params.outdir}/${meta.gca}` |
| versions_file | path | Versions YAML file tracking the script and Python versions |

## Parameters

### Required
- `params.update_registry`: Boolean flag to enable/disable process execution
- `params.asm_metadata`: JSON string (or an object in a `-params-file`) with [pymysql connection arguments](https://pymysql.readthedocs.io/en/latest/modules/connections.html) for the assembly metadata DB
- `params.outdir`: Output directory for results

Example on the command line:

```bash
--asm_metadata '{"host": "mysql-host", "port": 4527, "user": "rw_user", "password": "secret", "database": "gb_assembly_metadata"}'
```

Or in a `-params-file` (YAML):

```yaml
asm_metadata:
  host: mysql-host
  port: 4527
  user: rw_user
  password: secret
  database: gb_assembly_metadata
```

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process runs `bin/busco_metakeys_to_asmdb.py`, which:
1. Parses the BUSCO short summary file (modes `genome`, `euk_genome_met`, `euk_genome_min`)
2. Looks up `assembly_id` in the `assembly` table from `gca_chain` and `gca_version`
3. Inserts the metrics into `assembly_metrics` with `ON DUPLICATE KEY UPDATE`, so reruns overwrite existing values
4. Writes the parsed metrics to `<gca>_busco_genome_asmdb.json`

Metrics written: `assembly.busco`, `assembly.busco_version`, `assembly.busco_dataset`, `assembly.busco_mode`, `assembly.busco_completeness`, `assembly.busco_single_copy`, `assembly.busco_duplicated`, `assembly.busco_fragmented`, `assembly.busco_missing`, `assembly.busco_total` and `assembly.busco_erroneus`. `assembly.busco_erroneus` is only written when the BUSCO score line reports an `E:` value; otherwise it is omitted.

## Notes

- The assembly must already exist in the assembly metadata DB; the process fails otherwise
- Protein mode results are ignored by this process
- `asm_metadata` contains credentials and ends up in the task's `.command.sh` in the work directory; keep params files with it outside the repository
