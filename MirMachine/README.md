# ensembl-miRNA-nf
Automated MirMachine Based miRNA Detection workflow developed for genebuild team at Ensembl

## Deployment

### General

Default inputs to this Nextflow pipeline are specified in the params.config file. However, each can be overwritten when calling Nextflow:
Eg.
```
nextflow run main.nf --input <path to csv>
```
