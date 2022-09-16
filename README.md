# ensembl-genes-nf
Ensembl Genebuild NextFlow pipelines


Edit the config file specifying the scratch and the workDir


To run the pipeline in the cluster a csv with the list of db is needed. Moreover, in busco_pipeline.nf the db connections are specified but if the parameters are left blank, the db connections can be specified in the command. All the dbs in csv are supposed to have the same db connections.

BUSCO in genome mode
bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/busco_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv --genome_file genome.fa  --mode genome -w ../../work "
 
BUSCO in protein mode
bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/busco_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv --mode protein -w ../../work "
