# Busco Nextflow pipeline

Busco is a measure of completeness of genome assembly and annotation of the gene set. See the documentation for further details [Busco userguide](https://busco.ezlab.org/busco_userguide.html)

## Requirements

### Busco
We are using the Docker image available in https://hub.docker.com/r/ezlabgva/busco

### Ensembl dependencies
These are the Ensembl repositories required by this pipeline:

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

#### `--mode`
Select Busco mode, i.e. genome mode (assess a genome assembly), protein mode (assess a gene set) or both. By default, run both modes.

#### `--host`
The host name for the databases

#### `--port`
The port number of the host

#### `--user`
The read only username for the host. The password is expected to be empty.

#### `--enscode`
Path to the root directory containing the Perl repositories

#### `--outDir`
Path to the directory where to store the results of the pipeline

### Optional arguments

#### `--bioperl`
Path to the directory containing the BioPerl 1.6.924 library. If not provided, the value passed to `--enscode` will be used as root, i.e. `<enscode>/bioperl-1.6.924`.

#### `--cacheDir`
Path to the directory to use as cache for the intermediate files. If not provided, the value passed to `--outDir` will be used as root, i.e. `<outDir>/cache`.

#### `--project`
Output directory structure. Options: `ensembl`, `brc`. By default, `ensembl` project is selected.

#### `--files_latency`
Sleep time (in seconds) after the genome and proteins have been fetched. Needed by several file systems due to their internal latency. By default, 60 seconds.

### Pipeline configuration

#### Using the provided nextflow.config
We are using profiles to be able to run the pipeline on different HPC clusters. The default is `standard`.

* `standard`: uses LSF to run the compute heavy jobs. It expects the usage of `scratch` to use a low latency filesystem.
* `slurm`: uses SLURM to run the compute heavy jobs. It expects the usage of `scratch` to use a low latency filesystem.
* `local`: TODO


#### Using a local configuration file
You can use a local config with `-c` to finely configure your pipeline. All parameters can be configured, we recommend setting these ones as well:

* `process.scratch`: The patch to the scratch directory to use
* `workDir`: The directory where nextflow stores any file

### Running the different Busco modes
The default option is to run Busco in both genome and protein mode as follows:

#### Running BUSCO in LSF in genome mode
```bash
nextflow -C $ENSCODE/ensembl-genes-nf/nextflow.config run $ENSCODE/ensembl-genes-nf/workflows/busco_pipeline.nf --host <mysql_host> --port <mysql_port> --user <mysql_user> --enscode $ENSCODE --csvFile <csv_file_path> --mode genome

``` 
#### Running BUSCO in Slurm in protein mode
```bash
nextflow -C $ENSCODE/ensembl-genes-nf/nextflow.config run $ENSCODE/ensembl-genes-nf/workflows/busco_pipeline.nf --profile slurm --host <mysql_host> --port <mysql_port> --user <mysql_user> --enscode $ENSCODE --csvFile <csv_file_path> --mode protein
```

### Information about all the parameters

```bash
nextflow run ./ensembl-genes-nf/workflows/busco_pipeline.nf --help
```
