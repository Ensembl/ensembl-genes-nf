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

It is recommended that all the repositories are cloned into the same folder.

Remember that, following the instructions in [Ensembl's Perl API installation](http://www.ensembl.org/info/docs/api/api_installation.html), you will also need to have BioPerl v1.6.924 available in your system. If you do not, you can install it executing the following commands:

```bash
wget https://github.com/bioperl/bioperl-live/archive/release-1-6-924.zip
unzip release-1-6-924.zip
mv bioperl-live-release-1-6-924 bioperl-1.6.924
```

It is recommended to install it in the same folder as the Ensembl repositories.

## Running the pipeline


### Mandatory arguments

#### `--csvFile`
A file containing the list of databases you want to run Busco on. The databases need to have DNA.

#### `--host`
The host name for the databases

#### `--port`
The port number of the host

#### `--user`
The read only username for the host. The password is expected to be empty.

#### `--enscode`
Path to the root directory containing the Perl repositories

### Optional arguments

#### `--bioperl`
Path to the directory containing the BioPerl 1.6.924 library. If not provided, the value passed to `--enscode` will be used as root, i.e. `<enscode>/bioperl-1.6.924`.

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

## Running the pipeline with the Slurm profile

```bash
nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/workflows/omark_pipeline.nf -profile slurm --enscode $ENSCODE --csvFile dbname.csv -w ../../work
```
### Information about all the parameters

```bash
nextflow run ./ensembl-genes-nf/workflows/omark_pipeline.nf --help
```
