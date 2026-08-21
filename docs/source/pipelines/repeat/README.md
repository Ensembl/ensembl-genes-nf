# Repeat Pipeline

The **Repeat Pipeline** performs comprehensive repeat annotation of genome assemblies using a combination of de novo repeat detection and curated repeat libraries. Depending on the available resources and user configuration, the pipeline can generate species-specific RepeatModeler libraries, download existing libraries, annotate repeats with RepeatMasker, and perform complementary repeat detection using RED, DUST and Tandem Repeat Finder (TRF).

The pipeline is designed to support the Ensembl Genes annotation workflow by producing repeat annotations suitable for downstream analysis and FTP publication.

## Quick Links

- **[Input Specification](input.md)**
- **[Output Reference](output.md)**
- **[Parameter Reference](parameters.md)**
- **[Workflow Documentation](workflows/index.md)**
- **[Module Documentation](modules/index.md)**
- **[Troubleshooting](troubleshooting.md)**

---

## Pipeline Overview

The Repeat pipeline supports several repeat annotation strategies:

* Download genome assemblies from NCBI.
* Reuse existing RepeatModeler libraries when available.
* Generate new RepeatModeler libraries for assemblies without existing repeat models.
* Annotate repeats using RepeatMasker.
* Detect repeats independently using RED.
* Detect low-complexity regions using DUST.
* Detect tandem repeats using TRF.
* Upload generated repeat libraries and repeat annotations to the Ensembl FTP site.
* Record software versions for reproducibility.

Individual stages are enabled or disabled through pipeline parameters.

---

## Pipeline Workflow

The pipeline executes the following logical stages.

```text
Input CSV
    │
    ▼
FETCH_GENOME
    │
    ▼
FETCH_REPEAT_MODEL
    │
    ├──────────── Existing library ──────────────┐
    │                                            │
    ▼                                            ▼
GENERATE_REPEATMODELER_LIBRARY     CHECK_AND_DOWNLOAD_RMLIBRARY
    │                                            │
    └──────────────────────┬─────────────────────┘
                           ▼
                  RepeatModeler library
                           │
                           ▼
                    RUN_REPEATMASKER
                           │
        ┌──────────────────┼──────────────────┐
        ▼                  ▼                  ▼
      RUN_RED          RUN_DUST           RUN_TRF
        │                  │                  │
        └──────────────────┴──────────────────┘
                           ▼
              UPLOAD_REPEATS_INTO_FTP
                           │
                           ▼
           COLLECT_SOFTWARE_VERSIONS
```

---

## Main Components

### Genome Retrieval

The pipeline begins by downloading or locating the genome assembly specified in the input CSV.

**Module**

- [FETCH_GENOME](modules/fetch-genome.md)

---

### Repeat Library Preparation

When RepeatModeler library generation is enabled (`params.generate_lib`), the pipeline checks whether a species-specific repeat library already exists.

If one is available it is downloaded from the Ensembl FTP site.

Otherwise a new RepeatModeler library is generated and uploaded for future reuse.

**Modules**

- [FETCH_REPEAT_MODEL](modules/fetch-repeat-model.md)
- [GENERATE_REPEATMODELER_LIBRARY](modules/generate-repeatmodeler-library.md)
- [CHECK_AND_DOWNLOAD_RMLIBRARY](modules/check-and-download-rmlibrary.md)
- [UPLOAD_INTO_FTP](modules/upload-into-ftp.md)

---

### Repeat Annotation

RepeatMasker can be executed using the RepeatModeler library.

Additional repeat annotation methods can be enabled independently.

**Modules**

- [RUN_REPEATMASKER](modules/run-repeatmasker.md)
- [RUN_RED](modules/run-red.md)
- [RUN_DUST](modules/run-dust.md)
- [RUN_TRF](modules/run-trf.md)

Each method produces complementary repeat annotations.

---

### Result Publication

Generated repeat annotation files can optionally be uploaded to the Ensembl FTP site.

**Module**

- [UPLOAD_REPEATS_INTO_FTP] - (modules/upload-repeats-into-ftp.md)

---

### Reproducibility

All modules report software versions that are merged into a single `versions.yml` file.

**Module**

- [COLLECT_SOFTWARE_VERSIONS](modules/collect-software-versions.md)

---

## Execution Modes

The pipeline supports multiple execution modes depending on the supplied parameters.

### Complete Repeat Annotation

* Generate or retrieve RepeatModeler libraries.
* Run RepeatMasker.
* Run RED.
* Run DUST.
* Run TRF.
* Publish repeat annotations.

### RepeatMasker Only

Generate or retrieve RepeatModeler libraries and execute RepeatMasker only.

### RED Only

Execute RED directly on the genome assembly.

### DUST Only

Run DUST low-complexity masking.

### TRF Only

Detect tandem repeats using TRF.

Any combination of these analyses can be enabled.

---

## Pipeline Parameters

The generated parameter documentation provides the complete list of available configuration options.

See:

* `parameters.md`

---

## Input

The pipeline expects a CSV file describing the assemblies to process.

Typical columns include:

* `gca`
* `species_name`
* `genome_file`
* `repeatmasker_library`

See **Input Specification** for the complete format.

---

## Outputs

Depending on the enabled analyses, the pipeline generates:

* downloaded genome assemblies
* RepeatModeler libraries
* RepeatMasker annotations
* RED annotations
* DUST annotations
* TRF annotations
* FTP-ready output files
* software version reports

A detailed description of all outputs is available in **Output Reference**.

---

## Module Documentation

Detailed documentation is available for every module, including:

- Overview
- Inputs
- Outputs
- Parameters
- Implementation
- Dependencies
- Source

See the [Module Documentation](modules/index.md).

---

## Workflow Documentation

Workflow-level documentation describes how the pipeline orchestrates individual modules and how data flows between them.

See the [Workflow Documentation](workflows/index.md).

---

## Reproducibility

Every execution records software versions using the dedicated version collection module.

This produces a consolidated `versions.yml` file suitable for provenance tracking and reproducible analyses.

---

## Related Documentation

* Input Specification
* Output Reference
* Parameter Reference
* Module Documentation
* Workflow Documentation
* Troubleshooting
