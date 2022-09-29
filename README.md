# ensembl-genes-nf
Ensembl Genebuild NextFlow pipelines


Edit the config file specifying the scratch and the workDir


To run the pipeline in the cluster a csv with the list of db is needed. Moreover, in omark_pipeline.nf the db connections are specified but if the parameters are left blank, the db connections can be specified in the command. All the dbs in csv are supposed to have the same db connections.

bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/omark_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv  -w ../work"
