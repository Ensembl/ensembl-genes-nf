# Busco Nextflow pipeline

Busco is a easure of completeness of genome assembly and annotation of the gene set. See the documentation for further details https://busco.ezlab.org/busco_userguide.html
###Requirements

##Busco
We are using the Docker image available in https://hub.docker.com/r/ezlabgva/busco

### Perl EnsEMBL repositories you need to have

We recommend that you clone all the repositories into one directory
| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl | default | https://github.com/Ensembl/ensembl.git |
| ensembl-io | default | https://github.com/Ensembl/ensembl-io.git |

##Initialise the pipeline

Edit the config file specifying the scratch and the workDir


To run the pipeline in the cluster a csv with the list of database names is needed. Moreover, in busco_pipeline.nf the database connections are specified but they can be overwritten specifying the db connections in the command. All the dbs in csv are supposed to have the same db connections.

The default option is to run busco in both genome and protein mode

####BUSCO in genome mode
bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/busco_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv --genome_file genome.fa  --mode genome -w ../../work "
 
####BUSCO in protein mode
bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/busco_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv --mode protein -w ../../work "


Further information in the help
/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/busco_pipeline.nf --help
