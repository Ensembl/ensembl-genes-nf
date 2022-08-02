#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.db = ''
params.host = 'mysql-ens-sta-5'
params.port = '4684'
params.user = 'ensro'
params.pass = ''
params.work_dir = ''
params.enscode = ''
params.perl5lib = ''
params.busco_set = ''
params.busco_version = 'v5.3.2_cv1'
params.protein_file = 'translations.fa'
params.dump_params = ''
params.csvFile = ''
params.outDir = "/nfs/production/flicek/ensembl/genebuild/ftricomi/nextflow"
params.meta_file = ""
params.modules_path='/hps/software/users/ensembl/repositories/ftricomi/ensembl-genes-nf/modules.nf'
//Modules

include { queryDb } from params.modules_path
include { lowerLetters } from params.modules_path
include { getSpecies } from params.modules_path
include { getGca } from params.modules_path


//csvData = Channel.fromPath("${params.csvFile}").splitCsv(header: ['db'])
//csvData.into{csv_data1; csv_data2; csv_data3}

meta = Channel.fromPath("$params.meta_file}").collect()

/* dump canonical translations */
process fetchProteins {
  cpus 1
  memory { 6.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  clusterOptions = '-R "select[mem>2000] rusage[mem=2000]" -M2000'

  input:
  //tuple  val(host), val(port), val(user), val(db), val(dnahost), val(dnaport), val(dnauser), val(dnadb) from csv_data2
  
  val species_dir 
  val db

  storeDir "${params.outDir}/busco_score_test/${species_dir.trim()}/fasta/"

  output:
  path "translations.fa", emit: fasta
  val "${species_dir.trim()}", emit: output_dir
  val db, emit:dbname
  //beforeScript "source ${params.perl5lib}"
  //beforeScript "export PERL5LIB=${params.perl5lib}"
  script:
  //perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname $db -user ${params.user} -dnadbhost $dnahost -dnadbport $dnaport -dnadbname $dnadb -dnadbuser $dnauser -canonical_only 1 -file translations.fa  ${params.dump_params}
  """ 
  perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname $db -user ${params.user} -dnadbhost ${params.host} -dnadbport ${params.port} -dnadbname $db -dnadbuser ${params.user} -canonical_only 1 -file translations.fa  ${params.dump_params}
  """
}


/* run Busco in protein mode */
process runBusco {

  cpus 20
  memory { 40.GB * task.attempt }

  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  module 'singularity-3.7.0-gcc-9.3.0-dp5ffrp'
  container "ezlabgva/busco:${params.busco_version}"
  containerOptions "-B ${params.outDir}:/busco_wd"
  //runOptions = '--pull=always'
  clusterOptions = '-n 20 -R "select[mem>40000] rusage[mem=40000]" -M40000'

  input:
  file translations 
  //tuple val(host), val(port), val(user), val(db), val(dnahost), val(dnaport), val(dnauser), val(dnadb) from csv_data3
  val outdir 
  val db
 

  output:
  path "statistics/*.txt", emit: summary_file
  val outdir, emit:species_outdir

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.outDir}/busco_score_test/${outdir}/",  mode: 'copy' 

  script:
  println "${params.outDir}/busco_score_test/${outdir}/statistics/"
  
  //busco -f -i ${translations} --out busco_score_output --mode proteins -l ${params.busco_set} -c ${task.cpus}
  //singularity exec  --bind ${workflow.workDir}:/busco_wd /hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/busco-v5.1.2_cv1.simg  busco -f -i ${translations}  --mode proteins -l ${params.busco_set} -c ${task.cpus} -o statistics

  """
  busco -f -i ${translations}  --mode proteins -l ${params.busco_set} -c ${task.cpus} -o statistics --offline --download_path /nfs/production/flicek/ensembl/genebuild/ftricomi/busco_ftp/busco-data.ezlab.org/v5/data
  
  """
}

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

	queryDb(csvData.flatten(), meta)
	fetchProteins (queryDb.out.species_dir,queryDb.out.dbname)
	runBusco (fetchProteins.out.fasta.flatten(), fetchProteins.out.output_dir ,fetchProteins.out.dbname)
	lowerLetters (runBusco.out.species_outdir)
        getSpecies (lowerLetters.out.lower_name, lowerLetters.out.species_outdir)
        getGca (getSpecies.out.get_species, getSpecies.out.species_outdir, getSpecies.out.species)
        renameOutput (getGca.out.production_name, getGca.out.get_gca, getGca.out.species_outdir)
}
