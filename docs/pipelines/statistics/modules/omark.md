# OMark

## Overview

The `OMARK` process performs quality assessment of protein annotations using OMark (Orthology-based Marker assessment), analyzing OMAmer orthology assignments to evaluate annotation completeness and consistency.

## Process Details

- **Label**: `omamer`
- **Tag**: Uses genome assembly accession (`meta.gca`)
- **Publish Directory**: `${params.outdir}/${meta.gca}`
- **Max Forks**: 15 (limits parallel execution)

## Inputs

| Name | Type | Description |
|------|------|-------------|
| meta | val | Metadata map containing genome assembly info |
| omamer_file | path | OMAmer output file containing orthology assignments |

### Required Metadata Fields

- `gca`: Genome assembly accession

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| omark_output | tuple | Metadata, OMark text results, and all output files from `omark_output/` directory |
| versions_file | path | versions.yml file tracking OMark version |

## Parameters

### Required
- `params.outdir`: Output directory for results
- `params.omamer_database`: Path to OMAmer HOG database (used for reference during assessment)

### Optional
- `params.files_latency`: Delay after script execution (default handled by afterScript)

## Script Details

The process:
1. Executes `omark` on the OMAmer orthology file
2. Uses the OMAmer database as reference for quality assessment
3. Generates quality metrics including:
   - Completeness scores
   - Consistency metrics
   - Contamination detection
   - Fragment analysis
4. Outputs results to `omark_output/` directory
5. Captures OMark version from pip package information

## Dependencies

- OMark (Python package)
- OMAmer HOG database
- Python 3 with pip

## Output Files

The `omark_output/` directory contains:
- `*.txt` files: Summary statistics and quality metrics
- Additional output files: Detailed results and visualizations

## Quality Metrics

OMark provides several key quality indicators:
- **Completeness**: Percentage of expected orthologous groups present
- **Consistency**: Agreement between annotation and orthology
- **Contamination**: Detection of potential non-target sequences
- **Fragmentation**: Assessment of incomplete gene models

## Notes

- Results are published to a genome-specific subdirectory
- Maximum of 15 concurrent processes to manage computational resources
- The process includes a configurable sleep delay after completion to handle file system latency
- OMark version is extracted from pip package metadata
- Useful for validating gene annotation quality before database releas