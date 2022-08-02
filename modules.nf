/* Get Busco dataset using NSCBI taxonomy in meta table */
process BUSCODATASET {
  cpus 1
  memory { 2.GB * task.attempt }
  clusterOptions = '-R "select[mem>2000] rusage[mem=2000]" -M2000'
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  //beforeScript "source ${params.perl5lib}"
  //beforeScript "export PERL5LIB=${params.perl5lib}"
  beforeScript "export ${params.enscode}"
  beforeScript "source $ENSCODE/ensembl-genes-nf/supplementary_files/perl5lib.sh"
  input:
  tuple val(db)
  //file script

  output:
  val(db), emit:dbname
  stdout  emit:busco_dataset
  script:
  // get <Production name>/GCA
  """
  bash ${params.get_dataset_query} ${params.user}  ${params.host} ${params.port} $db
  """
  //mysql -N -u ${params.user}  -h ${params.host} -P ${params.port} -D $db < "${params.meta_file}"
  
}
/* Get species name and accession from meta table to build the output directory tree */
process SPECIESOUTDIR {
  cpus 1
  memory { 2.GB * task.attempt }
  clusterOptions = '-R "select[mem>2000] rusage[mem=2000]" -M2000'
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  tuple val(db)
  //file meta 
  val busco_dataset

  output:
  val(db), emit:dbname
  val busco_dataset, emit:busco_dataset
  stdout  emit:species_dir   

  script:
  // get <Production name>/GCA
  """
  mysql -N -u ${params.user}  -h ${params.host} -P ${params.port} -D $db < "${params.meta_query_file}"
  """
}

/* copy (and unzip) unmasked genome file */
process FETCHGENOME {
  cpus 1
  memory { 6.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  
  input:
  val species_dir 
  val db 
  val busco_dataset

  storeDir "${params.outDir}/${species_dir.trim()}/genome/"

  output:
  file "genome.fa", emit:fasta
  val "${species_dir.trim()}", emit:output_dir
  val db, emit:db_name
  val busco_dataset, emit:busco_dataset

  //check that the genome file is available 
  when:
  //file("/nfs/ftp/ensemblftp/ensembl/PUBLIC/pub/rapid-release/species/${species_dir.trim()}/genome").isDirectory()
  file("${params.genome_file}").isFile()

  script:
  """
  mkdir -p ${params.outDir}/${species_dir.trim()}/genome/
  cp "${params.genome_file}" ${params.outDir}/${species_dir.trim()}/genome/genome.fa
  """
  //cp /nfs/ftp/ensemblftp/ensembl/PUBLIC/pub/rapid-release/species/${species_dir.trim()}/genome/*-unmasked.fa.gz ${params.outDir}/busco_score_RR_NEW/${species_dir.trim()}/genome/genome.fa.gz
  //gzip -d -f ${params.outDir}/busco_score_RR_NEW/${species_dir.trim()}/genome/genome.fa.gz
  
}


process BUSCOGENOME {

  cpus 20
  memory { 60.GB * task.attempt }

  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  module 'singularity-3.7.0-gcc-9.3.0-dp5ffrp'
  container "ezlabgva/busco:${params.busco_version}"
  containerOptions "-B ${params.outDir}:/busco_wd"
  clusterOptions = '-n 20 -R "select[mem>60000] rusage[mem=60000]" -M60000'

  input:
  file genome 
  val output_dir 
  val db 
  val busco_dataset

  output:
  path "statistics/*.txt", emit:summary_file
  val output_dir, emit:species_outdir

  // ourdir is Production_name/GCA 
  publishDir "${params.outDir}/${output_dir}/",  mode: 'copy'

  script:
  println "${params.outDir}/${output_dir}/statistics/"

  """
  busco -f -i ${genome}  --mode genome -l ${busco_dataset}  -c ${task.cpus} -o statistics --offline --download_path /nfs/production/flicek/ensembl/genebuild/ftricomi/busco_ftp/busco-data.ezlab.org/v5/data
  """
}
/* Dump canonical translations */
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



/*ftp directory is Salmo_trutta/GCA_901001165.1/statistics/salmo_trutta_gca901001165v1_busco_short_summary.txt
in the following processes, summary file name is changed in <production name>_gca_busco_short_summary.txt

*/
process lowerLetters {
    // in: <Production name>/GCA
    // out: <production name>/gca
    input:
    val production_name
    output:
    stdout   emit:lower_name
    val production_name, emit:species_outdir

    """
    printf '$production_name' | tr '[A-Z]' '[a-z]' | tr . v
    """
}

process getSpecies {
    //in : <production name>/gcav1
    //out: <production name>
    input:
    val production_name
    val outdir

    output:
    stdout  emit:species
    val production_name, emit:get_species
    val outdir, emit:species_outdir
    """
    printf '$production_name' | cut -d'/' -f1
    """
}

process getGca {
    //in : <production name>/gcav1
    //out: gcav1
    input:
    val production_name 
    val outdir 
    val species_name 
    output:
    stdout emit:get_gca
    val outdir, emit:species_outdir
    val species_name, emit:production_name
    """
    printf '$production_name' | cut -d'/' -f2 | tr -d '_'
    """
}
