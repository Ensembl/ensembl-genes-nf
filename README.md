# Assembly Checker Nextflow pipeline
  
The pipeline aims to check the quality of an assembly downloading the genome sequences from the NCBI ftp and running Busco in genome mode. The assembly is defined good if the Busco score is above 90%.

### Requirements

### Busco
We are using the Docker image available in https://hub.docker.com/r/ezlabgva/busco

## Running the pipeline


### Mandatory options

#### csvFile
A file containing the list of "GCA,assembly name"

### Using the provided nextflow.config
We are using profiles to be able to run the pipeline on different HPC. The default is 'standard'

#### standard
Uses LSF to run the compute heavy jobs. It expects the usage of `scratch` to use a low latency filesystem

#### cluster
Uses SLURM to run the compute heavy jobs. It expects the usage of `scratch` to use a low latency filesystem

#### process.scratch
The patch to the scratch directory to use

#### workDir
The directory where nextflow stores any file

#### outDir
The directory to use to store the results of the pipeline


## Run the pipeline

```
nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/workflows/assembly_checker.nf --csvFile 
```

Further information in the help

```
nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/workflows/assembly_checker.nf --help
```
