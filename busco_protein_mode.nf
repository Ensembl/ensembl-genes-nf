#!/usr/bin/env nextflow
//nextflow.enable.dsl=2
params.db = ''
params.host = 'mysql-ens-sta-5'
params.port = '4684'
params.user = 'ensro'
params.pass = 'ensadmin'
params.work_dir = ''
params.enscode = ''
params.perl5lib = ''
params.busco_set = ''
params.busco_version = 'v5.3.2_cv1'
params.protein_file = 'translations.fa'
params.dump_params = ''
params.csvFile = ''

params.meta_file = "meta.sql"

//csvData = Channel.fromPath("${params.csvFile}").splitCsv(header: ['host', 'port', 'user', 'db', 'dnahost', 'dnaport', 'dnauser', 'dnadb'])

csvData = Channel.fromPath("${params.csvFile}").splitCsv(header: ['db'])
csvData.into{csv_data1; csv_data2; csv_data3}

meta = Channel.fromPath('meta.sql').collect()


process queryDb {
  cpus 1
  memory { 2.GB * task.attempt }
  //clusterOptions = "-R "select[mem>2000] rusage[mem=2000]" -M2000'
  // errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  //tuple val(host), val(port), val(user), val(db), val(dnahost), val(dnaport), val(dnauser), val(dnadb) from csv_data1
  tuple val(db) from csv_data1
  file meta from meta
  output:
  stdout  into query
  script:
  // get <Production name>/GCA
  """
  mysql -N -u ${params.user}  -h ${params.host} -P ${params.port} -D $db < "${meta}" 
  """
}


/* dump canonical translations */
process fetchProteins {
  cpus 1
  memory { 6.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  //tuple  val(host), val(port), val(user), val(db), val(dnahost), val(dnaport), val(dnauser), val(dnadb) from csv_data2
  tuple val(db) from csv_data2
  val species_dir from query

  storeDir "${workflow.outDir}/busco_score_RR/${species_dir.trim()}/fasta/"

  output:
  file "translations.fa" into fasta
  val "${species_dir.trim()}" into output_dir

  //beforeScript "source ${params.perl5lib}"
  //beforeScript "export PERL5LIB=${params.perl5lib}"
  script:
  println "VARIABILI"
  println "${species_dir.trim()}"
  //perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname $db -user ${params.user} -dnadbhost $dnahost -dnadbport $dnaport -dnadbname $dnadb -dnadbuser $dnauser -canonical_only 1 -file translations.fa  ${params.dump_params}
  """ 
  perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname $db -user ${params.user} -dnadbhost ${params.host} -dnadbport ${params.port} -dnadbname $db -dnadbuser ${params.user} -canonical_only 1 -file translations.fa  ${params.dump_params}
  """
}


/* run Busco in protein mode */
process runBusco {

  cpus 20
  memory { 20.GB * task.attempt }

  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  module 'singularity-3.7.0-gcc-9.3.0-dp5ffrp'
  container "ezlabgva/busco:${params.busco_version}"
  containerOptions "-B ${workflow.outDir}:/busco_wd"
  //runOptions = '--pull=always'

  input:
  file translations from fasta.flatten()
  //tuple val(host), val(port), val(user), val(db), val(dnahost), val(dnaport), val(dnauser), val(dnadb) from csv_data3
  tuple val(db) from csv_data3
  val outdir from output_dir

  output:
  path "statistics/*.txt" into summary_file
  val outdir into species_outdir
  // ourdir is Salmo_trutta (production name)
  publishDir "${workflow.outDir}/busco_score_RR/${outdir}/",  mode: 'copy' 

  script:
  println "${workflow.outDir}/busco_score_RR/${outdir}/statistics/"
  
  //busco -f -i ${translations} --out busco_score_output --mode proteins -l ${params.busco_set} -c ${task.cpus}
  //singularity exec  --bind ${workflow.workDir}:/busco_wd /hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/busco-v5.1.2_cv1.simg  busco -f -i ${translations}  --mode proteins -l ${params.busco_set} -c ${task.cpus} -o statistics

  """
  busco -f -i ${translations}  --mode proteins -l ${params.busco_set} -c ${task.cpus} -o statistics
  
  """
}

/*ftp directory is Salmo_trutta/GCA_901001165.1/statistics/salmo_trutta_gca901001165v1_busco_short_summary.txt 
in the following processes, summary file name is changed in <production name>_gca_busco_short_summary.txt

*/
process lowerLetters {
    // in: <Production name>/GCA
    // out: <production name>/gca
    input: 
    val production_name from species_outdir
    output: 
    stdout  into lower_case
    val production_name into species_outdir1  

    """ 
    printf '$production_name' | tr '[A-Z]' '[a-z]'
    """ 
} 

process modifyVersion {
    // in: <production name>/gca.1
    // out: <production name>/gcav1
    input:
    val production_name  from lower_case   
    val outdir from species_outdir1

    output:
    stdout into gca_version
    val outdir into species_outdir2

    """
    printf '$production_name' | tr . v
    """
}

process getSpecies {
    //in : <production name>/gcav1
    //out: <production name>
    input:
    val production_name  from gca_version
    val outdir from species_outdir2

    output:
    stdout into species
    val production_name  into get_species
    val outdir into species_outdir3
    """
    printf '$production_name' | cut -d'/' -f1
    """
}

process getGca {
    //in : <production name>/gcav1
    //out: gcav1
    input:
    val production_name  from get_species
    val outdir from species_outdir3
    output:
    stdout into get_gca
    val outdir into species_outdir4

    """
    printf '$production_name' | cut -d'/' -f2 | tr -d '_'
    """
}

process renameOutput {
    /*
	rename busco summary file in <production name>_gca_busco_short_summary.txt
	
    */
    input:
    val production_name from species
    val gca from get_gca
    val outdir from species_outdir4

    publishDir "${workflow.outDir}/busco_score_RR/${outdir}/",  mode: 'copy'
    output:
    //stdout into pippo
    //path busco into busco_path

    """
    mv -f ${workflow.outDir}/busco_score_RR/${outdir}/statistics/short_summary*  ${workflow.outDir}/busco_score_RR/${outdir}/statistics/${production_name.trim()}_${gca.trim()}_busco_short_summary.txt
    sed  -i '/genebuild/d' ${workflow.outDir}/busco_score_RR/${outdir}/statistics/${production_name.trim()}_${gca.trim()}_busco_short_summary.txt
    """
    //publishDir "${workflow.workDir}/busco_score_RR/${outdir}/*.txt", mode: "copy", pattern: 'busco_score_output/${x//_/[1]}_${x//_/[2]}_short_summary.txt'
}
