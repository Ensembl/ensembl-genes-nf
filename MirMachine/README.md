# ensembl-miRNA-nf
Automated MirMachine Based miRNA Detection workflow developed for genebuild team at Ensembl

## Deployment
The only requirements for running this workflow from a software perspective are [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Singularity/Apptainer](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages). 

Original development and testing was on EMBL-EBIs Codon HPC with further testing on a standalone Ubuntu 20.04 machine. Tested with Nextflow version 24.10.4 and Apptainer 1.3.4

### General

Default inputs to this Nextflow pipeline are specified in the params.config file. However, each can be overwritten when calling Nextflow:
Eg.
```
nextflow run main.nf \
  --input <path to tsv> \ # Tab separated file with columns titled 'Scientific Name' and 'Accession'. These must match what is in Ensembl FTP
  --outdir <path to output directory>
  --fasta_dir <path to where fastas are stored> # Used to speed up processing when fastas are on local file system
```
