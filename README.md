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

## Initialise the pipeline

Edit the config file specifying the scratch and the workDir


To run the pipeline in the cluster a csv with the list of database names is needed. Moreover, in omark_pipeline.nf the database connections are specified but they can be overwritten specifying the db connections in the command. All the dbs in csv are supposed to have the same db connections.

The default option is to run busco in both genome and protein mode

## Run the pipeline

```
bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/omark_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv  -w ../work"
```

Further information in the help

```
/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/omark_pipeline.nf --help# ensembl-genes-nf
```
