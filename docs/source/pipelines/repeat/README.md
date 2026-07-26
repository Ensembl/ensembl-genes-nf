# Repeat Pipeline Modules Documentation

This directory contains comprehensive documentation for all modules in the Ensembl genes repeat annotation pipeline.

## Module Overview

The repeat pipeline consists of 11 modules organised into functional categories:

### Library Retrieval Modules

1. **[check-and-download-rmlibrary](modules/check-and-download-rmlibrary.md)** – Downloads or validates RepeatModeler repeat libraries for the target assembly.
2. **[check-and-download-dfam](modules/check-and-download-dfam.md)** – Downloads and prepares Dfam repeat libraries when required.

### Repeat Annotation Modules

3. **[run-repeatmasker](modules/run-repeatmasker.md)** – Identifies and masks repetitive elements using RepeatMasker.
4. **[run-red](modules/run-red.md)** – Detects repeats de novo using RED.
5. **[merge-repeat-libraries](modules/merge-repeat-libraries.md)** – Combines multiple repeat libraries into a unified resource when required.

### Statistics & Reporting Modules

6. **[repeat-statistics](modules/repeat-statistics.md)** – Computes repeat annotation summary statistics.
7. **[repeat-summary](modules/repeat-summary.md)** – Generates summary reports for repeat annotation results.

### Database Operations Modules

8. **[populate-repeat-db](modules/populate-repeat-db.md)** – Populates Ensembl databases with repeat annotation metadata.
9. **[repeat-db-metadata](modules/repeat-db-metadata.md)** – Updates repeat-related metadata within the Ensembl core database.

### Resource Management Modules

10. **[clean-repeat-cache](modules/clean-repeat-cache.md)** – Removes cached intermediate files and temporary resources.
11. **[versions](modules/versions.md)** – Collects software version information for reproducibility.

> **Note:** Module names should be updated to match the final generated documentation if they differ from the current implementation.

---

## Pipeline Flow

The typical execution flow of the repeat pipeline:

```text
1. Repeat Library Preparation
   ├─> CHECK_AND_DOWNLOAD_RMLIBRARY
   └─> CHECK_AND_DOWNLOAD_DFAM

2. Repeat Annotation
   ├─> RUN_REPEATMASKER
   └─> RUN_RED

3. Result Processing
   ├─> MERGE_REPEAT_LIBRARIES
   ├─> REPEAT_STATISTICS
   └─> REPEAT_SUMMARY

4. Database Population
   ├─> POPULATE_REPEAT_DB
   └─> REPEAT_DB_METADATA

5. Cleanup
   └─> CLEAN_REPEAT_CACHE
```

---

## Module Categories by Function

### Repeat Library Management

* **RepeatModeler Library**: Retrieves or validates species-specific repeat libraries.
* **Dfam Library**: Downloads curated repeat family databases.

### Repeat Annotation

* **RepeatMasker**: Identifies and classifies known repetitive elements.
* **RED**: Detects repetitive regions using de novo sequence analysis.
* **Library Merging**: Combines multiple repeat resources for downstream analysis.

### Statistics & Reporting

* **Repeat Statistics**: Computes repeat coverage and annotation metrics.
* **Repeat Summary**: Produces summary reports for downstream quality assessment.

### Database Management

* **Populate Repeat Database**: Inserts repeat annotation results into Ensembl databases.
* **Repeat Metadata**: Updates database metadata and version information.

### Resource Management

* **Cleaning**: Removes temporary files and cached intermediate data.

---

## Key Dependencies

### External Tools

* **RepeatMasker** – Repeat annotation.
* **RepeatModeler** – De novo repeat family generation.
* **Dfam** – Curated repeat family database.
* **RED** – Rapid de novo repeat detection.

### Ensembl Dependencies

* **Ensembl Perl API** – Database interaction.
* **Ensembl Python libraries** – Pipeline utilities and metadata generation.
* **Ensembl analysis scripts** – Database loading and reporting.

### Databases

* **Ensembl Core Database** – Stores repeat annotations.
* **Dfam Database** – Reference repeat families.
* **RepeatModeler Libraries** – Species-specific repeat libraries.

---

## Common Parameters

Most modules use the following common parameters.

### Database Connection

* `params.host` – Database host.
* `params.port` – Database port.
* `params.user` – Database username.
* `params.password` – Database password.

### Paths

* `params.outdir` – Output directory.
* `params.cacheDir` – Cache directory.
* `params.enscode` – Ensembl code checkout.

### Execution Control

* `params.files_latency` – File-system synchronisation delay.
* `maxForks` – Maximum parallel processes.

---

## Caching Strategy

Several modules cache downloaded resources and intermediate files to avoid unnecessary recomputation.

Typical cached resources include:

* RepeatModeler libraries
* Dfam libraries
* RepeatMasker intermediate outputs
* RED intermediate files

This reduces download time and improves reproducibility when processing multiple assemblies.

---

## Conditional Execution

Some modules execute only when specific parameters are enabled.

Typical examples include:

* Use of RepeatModeler libraries versus Dfam libraries.
* Database population steps.
* Cleanup of working directories.
* Optional statistics generation.

Refer to the individual module documentation for the exact controlling parameters.

---

## Output Structure

Results are typically organised by genome assembly accession.

```text
${params.outdir}/
└── ${meta.gca}/
    ├── repeatmasker/
    ├── red/
    ├── repeat_library/
    ├── statistics/
    ├── reports/
    └── versions.yml
```

---

## Metadata Requirements

Most modules expect metadata maps containing fields such as:

* `gca` – Genome assembly accession.
* `dbname` – Ensembl core database.
* `production_name` – Species production name.
* `species_id` – Species identifier.

Additional metadata may be required by individual modules.

---

## Documentation Format

Each module documentation includes:

* **Overview**
* **Process Details**
* **Inputs**
* **Outputs**
* **Parameters**
* **Implementation Summary**
* **Dependencies**
* **Source**

---

## Version Tracking

Every module generates a `versions.yml` file recording software versions used during execution, including:

* RepeatMasker
* RepeatModeler
* RED
* Perl
* Python
* Database client versions

These files support reproducibility and troubleshooting.

---

## For More Information

* Refer to the individual module documentation for implementation details.
* Consult the main repeat pipeline documentation for workflow orchestration.
* See the Ensembl documentation for repeat annotation integration.
* Refer to the RepeatMasker, RepeatModeler, Dfam and RED documentation for tool-specific guidance.
