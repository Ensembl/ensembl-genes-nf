#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// db connection
params.db = ''
//params.host = 'mysql-ens-sta-5'
params.host = 'mysql-ens-genebuild-prod-6.ebi.ac.uk'
//params.port = '4684'
params.port = '4532'
params.user = 'ensro'
params.pass = ''
//repos
params.enscode = ''
params.perl5lib = ''
params.modules_path='/hps/software/users/ensembl/repositories/ftricomi/ensembl-genes-nf/modules.nf'

params.csvFile = ''
params.meta_query_file = '$ENSCODE/ensembl-genes-nf/supplementary_files/meta.sql'
params.get_dataset_query = '$ENSCODE/ensembl-genes-nf/supplementary_files/get_busco_dataset.sh'
params.outDir = "/nfs/production/flicek/ensembl/genebuild/ftricomi/nextflow/busco_score_test"
params.work_dir = ''
params.genome_file = ''

// Busco params
params.busco_set = ''
params.mode = ''
params.busco_version = 'v5.3.2_cv1'
params.download_path = '/nfs/production/flicek/ensembl/genebuild/ftricomi/busco_ftp/busco-data.ezlab.org/v5/data'
params.dump_params = ''
params.meta_file = "$ENSCODE/ensembl-genes-nf/supplementary_files/meta.sql"


//Modules

include { BUSCODATASET } from params.modules_path
include { SPECIESOUTDIR } from params.modules_path
include { FETCHGENOME } from params.modules_path
include { FETCHPROTEINS } from params.modules_path
include { BUSCOGENOME } from params.modules_path
include { BUSCOPROTEIN } from params.modules_path
include { GETSPECIESNAME } from params.modules_path
include { GETGCA } from params.modules_path
include { OUTPUT } from params.modules_path
include { INPUT_CHECK } from params.modules_path

params.help = false

 // print usage
if (params.help) {
  log.info ''
  log.info 'Pipeline to run Busco score in protein and/or genome mode'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C ensembl-genes-nf/nextflow.config run ensembl-genes-nf/busco_pipeline.nf --enscode --csvFile --genome_file --mode'
  log.info ''
  log.info 'Options:'
  log.info '  --host                    Db host server '
  log.info '  --port                    Db port  '
  log.info '  --user                    Db user  '
  log.info '  --enscode                 Enscode path '
  log.info '  --perl5lib                PERL5LIB '
  log.info '  --outDir                  Output directory '
  log.info '  --csvfile                 Path for the csv containing the db name'
  log.info '  --mode                    Busco mode: genome or protein'
  log.info '  --genome_file      	FASTA genome file (unmasked)'
  log.info '  --cpus INT	        Number of CPUs to use. Default 1.'
  exit 1
}

if( !params.host) {
  exit 1, "Undefined --host parameter. Please provide the server host for the db connection"
}

if( !params.port) {
  exit 1, "Undefined --port parameter. Please provide the server port for the db connection"
}
if( !params.user) {
  exit 1, "Undefined --user parameter. Please provide the server user for the db connection"
}

if( !params.enscode) {
  exit 1, "Undefined --enscode parameter. Please provide the enscode path"
}
if( !params.perl5lib) {
  exit 1, "Undefined --perl5lib parameter. Please provide the PERL5LIB path"
}
if( !params.outDir) {
  exit 1, "Undefined --outDir parameter. Please provide the output directory's path"
}
if( !params.mode) {
  exit 1, "Undefined --mode parameter. Please define Busco running mode"
}
csvFile = file(params.csvfile)
if( !csvFile.exists() ) {
  exit 1, "The specified csv file does not exist: ${params.csvfile}"
}

workflow{
        csvData = Channel.fromPath("${params.csvFile}").splitCsv(header: ['db'])
        mode = Channel.from("${params.mode}").view()
        println "${mode}"
        BUSCODATASET (csvData.flatten())
	SPECIESOUTDIR (BUSCODATASET.out.dbname, BUSCODATASET.out.busco_dataset, params.mode)
        SPECIESOUTDIR.out.branch {
                        protein: it[3] == 'protein'
                        genome: it[3] == 'genome'
             }.set { ch_mode }
        //INPUT_CHECK ( SPECIESOUTDIR.out.species_dir, SPECIESOUTDIR.out.dbname, SPECIESOUTDIR.out.busco_dataset).branch { 
       // 		protein: params.mode == 'protein'
	//			println "CSSOOO"
	//			return [species_dir, dbname, busco_dataset]
        //		genome: params.mode == 'genome'        		            
         //    }.set { ch_mode }
        //println "$ch_mode"
        ch_mode.protein.view {"puzzaaa"} 
        
        FETCHGENOME (ch_mode.genome)
        BUSCOGENOME (FETCHGENOME.out.fasta.flatten(), FETCHGENOME.out.output_dir, FETCHGENOME.out.db_name, FETCHGENOME.out.busco_dataset)
        
        FETCHPROTEINS (ch_mode.protein)
        BUSCOPROTEIN (FETCHPROTEINS.out.fasta.flatten(), FETCHPROTEINS.out.output_dir, FETCHPROTEINS.out.db_name, FETCHPROTEINS.out.busco_dataset)
         
        //if(ch_mode.genome){
        //GETSPECIESNAME(BUSCOPROTEIN.out.species_outdir)
        //GETGCA(GETSPECIESNAME.out.species_outdir, GETSPECIESNAME.out.species_name)
        //OUTPUT(GETGCA.out.species_name, GETGCA.out.get_gca, GETGCA.out.species_outdir)
//}
}
