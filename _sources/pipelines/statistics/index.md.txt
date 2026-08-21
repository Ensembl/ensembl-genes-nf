# Statistics Pipeline

The **Statistics Pipeline** generates quality metrics and statistics for Ensembl gene annotations and assemblies.

## Overview

The pipeline provides:

- BUSCO completeness assessment
- OMArk annotation quality assessment
- Ensembl database statistics
- Metadata generation

## Documentation

```{toctree}
:maxdepth: 1

README
input
output
parameters
workflows/index
modules/index
troubleshooting
```

## Quick Start

```bash
nextflow run main.nf \
  --csvFile input.csv \
  --run_busco_core \
  --run_omark \
  --outdir results
```