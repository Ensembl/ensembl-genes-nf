# ensembl-genes-nf
Ensembl Genebuild NextFlow pipelines

To run the pipeline in the cluster 

BUSCO in genome mode
bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/busco_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv --genome_file genome.fa  --mode genome -w ../../work "
 
BUSCO in protein mode
bsub -Is "/hps/software/users/ensembl/genebuild/bin/nextflow -C ./ensembl-genes-nf/nextflow.config run ./ensembl-genes-nf/busco_pipeline.nf --enscode $ENSCODE --csvFile dbname.csv --mode protein -w ../../work "
