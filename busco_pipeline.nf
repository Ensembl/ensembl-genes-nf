#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// db connection
params.db = ''
params.host = 'mysql-ens-genebuild-prod-6.ebi.ac.uk'
params.port = '4532'
params.user = 'ensro'
params.pass = ''
//repos
params.enscode = ''
params.modules_path="${projectDir}/modules.nf"

params.csvFile = ''
params.meta_query_file = "${projectDir}/supplementary_files/meta.sql"
params.get_dataset_query = "${projectDir}/supplementary_files/get_busco_dataset.sh"
params.ortho_list = "${projectDir}/supplementary_files/busco_dataset.txt"
params.outDir = "${workDir}"
params.genome_file = ''
params.help = false

// Busco params
params.mode = ['protein', 'genome']
params.busco_version = 'v5.3.2_cv1'
params.download_path = '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/data/busco_data/data'

// Params for dumping the translations. Default is canonical transcripts only
// You should use '--nocanonical_only' if you want all proteins to be dumped
params.dump_params = '--canonical_only'


//Modules

include { BUSCODATASET } from params.modules_path
include { SPECIESOUTDIR } from params.modules_path
include { FETCHGENOME } from params.modules_path
include { FETCHPROTEINS } from params.modules_path
include { BUSCOGENOME } from params.modules_path
include { BUSCOPROTEIN } from params.modules_path
include { BUSCOGENOMEOUTPUT } from params.modules_path
include { BUSCOPROTEINOUTPUT } from params.modules_path

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
  log.info '  --host STR                   Db host server '
  log.info '  --port INT                   Db port  '
  log.info '  --user STR                   Db user  '
  log.info '  --enscode STR                Enscode path '
  log.info '  --outDir STR                 Output directory. Default is workDir'
  log.info '  --csvFile STR                Path for the csv containing the db name'
  log.info '  --mode STR                   Busco mode: genome or protein, default is to run both'
  log.info '  --genome_file STR            Unmasked FASTA genome file (genome mode)'
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
if( !params.mode) {
  exit 1, "Undefined --mode parameter. Please define Busco running mode"
}

csvFile = file(params.csvFile)
if( !csvFile.exists() ) {
  exit 1, "The specified csv file does not exist: ${params.csvfile}"
}
busco_mode = []
if (params.mode instanceof java.lang.String) {
  busco_mode = [params.mode]
}
else {
  busco_mode = params.mode
}
workflow{
        csvData = Channel.fromPath(params.csvFile).splitCsv()
        buscoModes = Channel.fromList(busco_mode)
        
        BUSCODATASET (csvData.flatten())
        SPECIESOUTDIR (BUSCODATASET.out.dbname, BUSCODATASET.out.busco_dataset)
        SPECIESOUTDIR.out.combine(buscoModes).branch {
                        protein: it[3] == 'protein'
                        genome: it[3] == 'genome'
             }.set { ch_mode }
        
        FETCHGENOME (ch_mode.genome)
        BUSCOGENOME (FETCHGENOME.out.fasta.flatten(), FETCHGENOME.out.output_dir, FETCHGENOME.out.db_name, FETCHGENOME.out.busco_dataset)
        BUSCOGENOMEOUTPUT(BUSCOGENOME.out.species_outdir)        
        FETCHPROTEINS (ch_mode.protein)
        BUSCOPROTEIN (FETCHPROTEINS.out.fasta.flatten(), FETCHPROTEINS.out.output_dir, FETCHPROTEINS.out.db_name, FETCHPROTEINS.out.busco_dataset)
        BUSCOPROTEINOUTPUT(BUSCOPROTEIN.out.species_outdir) 
}
