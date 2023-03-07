# OMArk Nextflow pipeline
  
OMArk is a software of proteome (protein-coding gene repertoire) quality assessment. It provides measure of proteome completeness, characterize all protein coding genes in the light of existing homologs, and identify the presence of contamination from other species.
Further information available in the official repo https://github.com/DessimozLab/OMArk

### Requirements

## OMArk
We built a  Docker image available in the cluster as singularity

### Perl EnsEMBL repositories you need to have

We recommend that you clone all the repositories into one directory
| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl | default | https://github.com/Ensembl/ensembl.git |
| ensembl-analysis | default | https://github.com/Ensembl/ensembl-analysis.git |
| ensembl-io | default | https://github.com/Ensembl/ensembl-io.git |


## Running the pipeline


### Mandatory options

#### csvFile
A file containing the list of databases you want to run Busco on. The databases need to have DNA.

#### host
The host name for the databases

#### port
The port number of the host

#### user
The read only username for the host. The password is expected to be empty.

#### enscode
The directory containing the Perl repositories


### Using the provided nextflow.config
We are using profiles to be able to run the pipeline on different HPC. The default is 'standard'

#### standard
Uses LSF to run the compute heavy jobs. It expects the usage of `scratch` to use a low latency filesystem

#### cluster
Uses SLURM to run the compute heavy jobs. It expects the usage of `scratch` to use a low latency filesystem


### Using a local config
You can use a local config with `-c` to finely configure your pipeline. All parameters can be configured, we recommend setting the ones mentionned below.

#### process.scratch
The patch to the scratch directory to use

#### workDir
The directory where nextflow stores any file

#### outDir
The directory to use to store the results of the pipeline

## Run the pipeline

```
nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/workflows/omark_pipeline.nf -profile slurm --enscode $ENSCODE --csvFile dbname.csv -w ../../work
```
### Information about all the parameters

```
nextflow run ./ensembl-genes-nf/workflows/omark_pipeline.nf --help
```
