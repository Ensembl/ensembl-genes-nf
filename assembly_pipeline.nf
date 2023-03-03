#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//repos
params.enscode = ''
params.modules_path='/hps/software/users/ensembl/repositories/ftricomi/ensembl-genes-nf/modules.nf'
params.meta_query_file = '$ENSCODE/ensembl-genes-nf/supplementary_files/meta.sql'
params.csvFile = ''
params.outDir = "/nfs/production/flicek/ensembl/genebuild/ftricomi/nextflow/busco_genome"
params.genome_file = ''

// Busco params
params.busco_version = 'v5.3.2_cv1'
params.download_path = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/data/busco_data/data'

//Modules

include { PROCESSASSEMBLY } from params.modules_path
include { BUSCOGENOME } from params.modules_path

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
csvFile = file(params.csvFile)
if( !csvFile.exists() ) {
  exit 1, "The specified csv file does not exist: ${params.csvfile}"
}

workflow{
        //csvData = Channel.fromPath("${params.csvFile}").splitCsv(sep:',').multiMap { it ->
        //gca: [it[0]]
        //assembly_name: [it[1]]
        //}
	    //PROCESSASSEMBLY (csvData.gca, csvData.assembly_name)
        csvData = Channel.fromPath(params.csvFile).splitCsv(sep:',')
        PROCESSASSEMBLY (csvData)
        BUSCOGENOME (PROCESSASSEMBLY.out.gca, PROCESSASSEMBLY.out.genome_file.flatten())
}
