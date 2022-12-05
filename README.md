# Busco Nextflow pipeline

Busco is a measure of completeness of genome assembly and annotation of the gene set. See the documentation for further details [Busco userguide](https://busco.ezlab.org/busco_userguide.html)

## Requirements

### Busco
We are using the Docker image available in https://hub.docker.com/r/ezlabgva/busco

### Perl EnsEMBL repositories you need to have
We recommend that you clone all the repositories into one directory
| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl | default | https://github.com/Ensembl/ensembl.git |
| ensembl-analysis | default | https://github.com/Ensembl/ensembl-analysis.git |


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


### Running the different Busco modes
The default option is to run busco in both genome and protein mode

#### BUSCO in genome mode

```
/hps/software/users/ensembl/genebuild/bin/nextflow run ./ensembl-genes-nf/busco_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv --genome_file genome.fa  --mode genome -w ../../work
``` 
#### BUSCO in protein mode

```
/hps/software/users/ensembl/genebuild/bin/nextflow run ./ensembl-genes-nf/busco_pipeline.nf -profile slurm --enscode $ENSCODE --csvFile dbname.csv --mode protein -w ../../work
```

### Information about all the parameters

```
/hps/software/users/ensembl/genebuild/bin/nextflow run ./ensembl-genes-nf/busco_pipeline.nf --help
```
