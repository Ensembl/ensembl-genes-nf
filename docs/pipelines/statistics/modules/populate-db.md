# Populate Database

## Overview

The `POPULATE_DB` process executes SQL files against Ensembl core databases to populate them with statistics, metadata, or other computed information.

## Process Details

- **Label**: `default`
- **Tag**: Uses database name (`meta.dbname`)
- **Conditional Execution**: Only runs when `params.apply_ensembl_stats` OR `params.apply_ensembl_beta_metakeys` is true

## Inputs

| Name | Type | Description |
|------|------|-------------|
| meta | val | Metadata map containing database connection details |
| sql_file | path | SQL file to execute against the database |

### Required Metadata Fields

- `dbname`: Target database name

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| versions_file | path | Optional versions.yml file tracking MySQL version |

## Parameters

### Required
- `params.mysql_ensadmin`: Path to MySQL admin script/command
- `params.host`: Database host
- `params.apply_ensembl_stats`: Boolean flag to enable statistics population
- `params.apply_ensembl_beta_metakeys`: Boolean flag to enable beta metakeys population

## Script Details

The process:
1. Executes the mysql_ensadmin script with host and database name
2. Redirects the SQL file as input to populate the database
3. Captures MySQL version information
4. Generates a versions file

### Command Structure
```bash
${params.mysql_ensadmin}/${params.host} ${meta.dbname} < ${sql_file}
```

## Dependencies

- MySQL client
- `mysql_ensadmin` script (Ensembl admin utility)
- Database write access with appropriate privileges

## Use Cases

This process is used to:
- Load gene and transcript statistics into core databases
- Insert computed metadata keys
- Apply database schema updates
- Populate summary tables

## Conditional Execution

The process only runs when at least one of these conditions is true:
- `params.apply_ensembl_stats == true`: Apply statistical summaries
- `params.apply_ensembl_beta_metakeys == true`: Apply beta metadata keys

## Notes

- The versions file is marked as optional
- SQL execution is performed via the mysql_ensadmin utility for standardized database administration
- The process assumes SQL files are pre-validated and safe to execute
- No rollback mechanism is provided; ensure SQL files are tested before production use