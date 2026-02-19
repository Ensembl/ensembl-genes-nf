# Ensembl Genes Nextflow Pipelines

Welcome to the documentation for the Ensembl Genes Nextflow pipelines! This repository contains production-grade bioinformatics workflows for genome annotation quality assessment and statistics generation.

## Overview

The Ensembl Genes pipelines provide automated workflows for:

- **Quality Assessment**: Evaluate genome assemblies and annotations using industry-standard metrics
- **Statistics Generation**: Compute comprehensive gene set statistics for Ensembl databases
- **Metadata Management**: Generate and apply metadata keys to Ensembl core databases

## Available Pipelines

### Statistics Pipeline

The statistics pipeline performs quality control and generates metrics for genome annotations, including:

- **BUSCO Analysis**: Assess completeness of genome assemblies and gene sets
- **OMArk Evaluation**: Measure proteome completeness based on conserved genes
- **Ensembl Statistics**: Generate core statistics for Ensembl databases

[Explore the Statistics Pipeline →](pipelines/statistics/overview.md){ .md-button .md-button--primary }

## Key Features

✅ **Reproducible**: Built with Nextflow for portable, reproducible workflows  
✅ **Scalable**: Designed to run on local machines, HPC clusters, or cloud platforms  
✅ **Comprehensive**: Multiple quality metrics and statistics in one pipeline  
✅ **Flexible**: Configurable parameters for different use cases  
✅ **Production-Ready**: Used by Ensembl for production genome annotations  

## Quick Links

- [Quick Start Guide](pipelines/statistics/quickstart.md)
- [Input Format](pipelines/statistics/input.md)
- [Parameters Reference](pipelines/statistics/parameters.md)
- [GitHub Repository](https://github.com/Ensembl/ensembl-genes-nf)

## Requirements

- **Nextflow**: Version 23.04 or later
- **Container Engine**: Docker, Singularity, or Podman
- **Java**: Version 11 or later (for Nextflow)

## Getting Started

```bash
# Clone the repository
git clone https://github.com/Ensembl/ensembl-genes-nf.git
cd ensembl-genes-nf/pipelines/statistics

# Run with example data
nextflow run main.nf -profile docker --csvFile data/example.csv --outdir results
```

## Support

For questions, issues, or contributions:

- **Issues**: [GitHub Issues](https://github.com/Ensembl/ensembl-genes-nf/issues)
- **Email**: [Ensembl Helpdesk](http://www.ensembl.org/Help/Contact)

## License

This project is licensed under the Apache License 2.0.
