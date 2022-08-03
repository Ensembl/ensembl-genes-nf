#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// db connection
params.db = ''
params.host = 'mysql-ens-sta-5'
params.port = '4684'
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
params.modules_path='/hps/software/users/ensembl/repositories/ftricomi/ensembl-genes-nf/modules.nf'


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




process renameOutput {
    /*
	rename busco summary file in <production name>_gca_busco_short_summary.txt
    */
    input:
    val production_name 
    val gca 
    val outdir 

    publishDir "${params.outDir}/busco_score_RR/${outdir}/",  mode: 'copy'

    """
    mv -f ${params.outDir}/busco_score_test/${outdir}/statistics/short_summary*  ${params.outDir}/busco_score_test/${outdir}/statistics/${production_name.trim()}_${gca.trim()}_busco_short_summary.txt
    sed  -i '/genebuild/d' ${params.outDir}/busco_score_test/${outdir}/statistics/${production_name.trim()}_${gca.trim()}_busco_short_summary.txt
    """
}
workflow {
        csvData = Channel.fromPath("${params.csvFile}").splitCsv(header: ['db'])
        BUSCODATASET (csvData.flatten())
	SPECIESOUTDIR (BUSCODATASET.out.dbname, BUSCODATASET.out.busco_dataset)
        if(params.mode == 'genome'){
          	FETCHGENOME (SPECIESOUTDIR.out.species_dir, SPECIESOUTDIR.out.dbname, SPECIESOUTDIR.out.busco_dataset)
		BUSCOGENOME (FETCHGENOME.out.fasta.flatten(), FETCHGENOME.out.output_dir, FETCHGENOME.out.db_name, FETCHGENOME.out.busco_dataset)
                GETSPECIESNAME(BUSCOGENOME.out.species_outdir)
        }
        if(params.mode == 'protein'){
                FETCHPROTEINS (SPECIESOUTDIR.out.species_dir, SPECIESOUTDIR.out.dbname, SPECIESOUTDIR.out.busco_dataset)
                BUSCOPROTEIN (FETCHPROTEINS.out.fasta.flatten(), FETCHPROTEINS.out.output_dir, FETCHPROTEINS.out.db_name, FETCHPROTEINS.out.busco_dataset)
                GETSPECIESNAME(BUSCOPROTEIN.out.species_outdir)
        }
        GETGCA(GETSPECIESNAME.out.species_outdir, GETSPECIESNAME.out.species_name)
        OUTPUT(GETGCA.out.species_name, GETGCA.out.get_gca, GETGCA.out.species_outdir)
}
