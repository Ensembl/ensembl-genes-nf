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
params.modules_path='/hps/software/users/ensembl/repositories/ftricomi/ensembl-genes-nf/modules.nf'
params.meta_query_file = '$ENSCODE/ensembl-genes-nf/supplementary_files/meta.sql'
params.csvFile = ''
params.outDir = "/nfs/production/flicek/ensembl/genebuild/ftricomi/nextflow/omark"
params.genome_file = ''

// Busco params
params.omamer_database = '/nfs/production/flicek/ensembl/genebuild/ftricomi/omark_database/OMAmerDB/LUCA_MinFamSize6_OR_MinFamComp05_A21_k6.h5 '
params.dump_params = ''


//Modules

include { SPECIESOUTDIR } from params.modules_path
include { FETCHPROTEINS } from params.modules_path
include { CREATEOMAMER } from params.modules_path
include { RUNOMARK } from params.modules_path

params.help = false

 // print usage
if (params.help) {
  log.info ''
  log.info 'Pipeline to run OMArk score measuring proteome (protein-coding gene repertoire) quality assessment'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C ensembl-genes-nf/nextflow.config run ensembl-genes-nf/omark_pipeline.nf --enscode --csvFile '
  log.info ''
  log.info 'Options:'
  log.info '  --host                    Db host server '
  log.info '  --port                    Db port  '
  log.info '  --user                    Db user  '
  log.info '  --enscode                 Enscode path '
  log.info '  --outDir                  Output directory '
  log.info '  --csvFile                 Path for the csv containing the db name'
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
if( !params.outDir) {
  exit 1, "Undefined --outDir parameter. Please provide the output directory's path"
}
csvFile = file(params.csvFile)
if( !csvFile.exists() ) {
  exit 1, "The specified csv file does not exist: ${params.csvfile}"
}

workflow{
        csvData = Channel.fromPath("${params.csvFile}").splitCsv()
	SPECIESOUTDIR (csvData.flatten())
        FETCHPROTEINS (SPECIESOUTDIR.out)
        CREATEOMAMER (FETCHPROTEINS.out.fasta.flatten(), FETCHPROTEINS.out.output_dir, FETCHPROTEINS.out.db_name)
        RUNOMARK (CREATEOMAMER.out.omamer_file, CREATEOMAMER.out.species_outdir)
         
}
