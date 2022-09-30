def get_species_name(name) {
  final m = name.tr('[A-Z]','[a-z]').tr('.','v').split('/')[0]
  return m
}

def get_gca(name) {
  final m = name.tr('[A-Z]', '[a-z]').tr('.', 'v').replaceAll("_","").split('/').getAt(1)
  return m
}
def get_gca_anno(name) {
  final m = name.split('/').getAt(1)
  return m
}

def concatString(string1, string2, string3){
 return string1 + '_'+string2 + '_'+string3
}

/* Get Busco dataset using NSCBI taxonomy in meta table */
process BUSCODATASET {
  cpus 1
  memory { 2.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  beforeScript "export ENSCODE=${params.enscode}"
  input:
  val db

  output:
  val db, emit:dbname
  stdout  emit:busco_dataset
  script:
  // get <Production name>/GCA
  """
  bash ${params.get_dataset_query} ${params.user} ${params.host} ${params.port} ${db} ${params.ortho_list}
  """
}
/* Get species name and accession from meta table to build the output directory tree */
process SPECIESOUTDIR {
  cpus 1
  memory { 2.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  val db
  val busco_dataset
  val mode
  output:
  tuple stdout, val(db), val(busco_dataset),val(mode)
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
  tuple val(species_dir),val(db), val(busco_dataset), val(mode)

  storeDir "${params.outDir}/${species_dir.trim()}/genome/"

  output:
  path "genome.fa", emit:fasta
  val "${species_dir.trim()}", emit:output_dir
  val db, emit:db_name
  val "${busco_dataset.trim()}", emit:busco_dataset
  //check that the genome file is available 
  //when:
  //file("/nfs/ftp/ensemblftp/ensembl/PUBLIC/pub/rapid-release/species/${species_dir.trim()}/genome").isDirectory()
  // file("${params.genome_file}").isFile()
  //file(/nfs/production/flicek/ensembl/genebuild/ftricomi/anno_annotation/get_gca("${species_dir.trim()}")).isDirectory()
  //a=get_gca("${species_dir.trim()}")
  script:
  """
  mkdir -p ${params.outDir}/${species_dir.trim()}/genome/
  cp /hps/nobackup/flicek/ensembl/genebuild/ftricomi/anno_annotations/from_nfs/${get_gca_anno("${species_dir.trim()}")}/*_reheadered_toplevel.fa ${params.outDir}/${species_dir.trim()}/genome/genome.fa
  """
  //cp "${params.genome_file}" ${params.outDir}/${species_dir.trim()}/genome/genome.fa
  
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
  containerOptions "-B ${params.outDir}/:/busco_wd"

  input:
  file genome
  val outdir
  val db
  val busco_dataset

  output:
  path "genome/*.txt", emit: summary_file
  val outdir, emit:species_outdir

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.outDir}/${outdir}/",  mode: 'copy'

  script:
  """
  busco -f -i ${genome} -o genome --mode genome -l ${busco_dataset} -c ${task.cpus} --offline --download_path ${params.download_path}
  """
}


 process BUSCOGENOMEOUTPUT {
     /*
         rename busco summary file in <production name>_gca_genome_busco_short_summary.txt
     */
     input:
     val outdir

     publishDir "${params.outDir}/${outdir}/",  mode: 'copy'
     script:
     """
     mkdir -p  ${params.outDir}/${outdir}/statistics
     sed  -i '/genebuild/d' ${params.outDir}/${outdir}/genome/short_summary* 
     mv -f ${params.outDir}/${outdir}/genome/short_summary* ${params.outDir}/${outdir}/statistics/${concatString(get_species_name("${outdir.trim()}"),get_gca("${outdir.trim()}"),'genome_busco_short_summary.txt')}
     """
 }



/* Dump canonical translations */
process FETCHPROTEINS {
  cpus 1
  memory { 6.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  tuple val(species_dir),val(db), val(busco_dataset), val(mode)
  storeDir "${params.outDir}/${species_dir.trim()}/fasta/"

  output:
  path "translations.fa", emit: fasta
  val "${species_dir.trim()}", emit: output_dir
  val db, emit:db_name
  val "${busco_dataset.trim()}", emit:busco_dataset  

  beforeScript "ENSCODE=${params.enscode} source ${projectDir}/supplementary_files/perl5lib.sh"
 
  script:
  """
  perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname $db -user ${params.user} -dnadbhost ${params.host} -dnadbport ${params.port} -dnadbname $db -dnadbuser ${params.user} -canonical_only 1 -file translations.fa  ${params.dump_params}
  """
}


/* run Busco in protein mode */
process BUSCOPROTEIN {

  cpus 20
  memory { 40.GB * task.attempt }

  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  module 'singularity-3.7.0-gcc-9.3.0-dp5ffrp'
  container "ezlabgva/busco:${params.busco_version}"
  containerOptions "-B ${params.outDir}:/busco_wd"

  input:
  file translations
  val outdir
  val db
  val busco_dataset

  output:
  path "fasta/*.txt", emit: summary_file
  val outdir, emit:species_outdir

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.outDir}/${outdir}/",  mode: 'copy'

  script:
  """
  busco -f -i ${translations}  --mode proteins -l ${busco_dataset} -c ${task.cpus} -o fasta --offline --download_path ${params.download_path}
  """
  //sed  -i '/genebuild/d' ${params.outDir}/${outdir}/fasta/short_summary*
  //mkdir -p ${params.outDir}/${outdir}/statistics 
  //mv -f ${params.outDir}/${outdir}/fasta/short_summary*  ${params.outDir}/${outdir}/statistics/${concatString(get_species_name("${outdir.trim()}"),get_gca("${outdir.trim()}"),'busco_short_summary.txt')}
  //"""
}

process BUSCOPROTEINOUTPUT {
     /*
         rename busco summary file in <production name>_gca_genome_busco_short_summary.txt
     */
     input:
     val outdir

     publishDir "${params.outDir}/${outdir}/",  mode: 'copy'
     script:
     """
     mkdir -p  ${params.outDir}/${outdir}/statistics
     sed  -i '/genebuild/d' ${params.outDir}/${outdir}/fasta/short_summary*
     mv -f ${params.outDir}/${outdir}/fasta/short_summary* ${params.outDir}/${outdir}/statistics/${concatString(get_species_name("${outdir.trim()}"),get_gca("${outdir.trim()}"),'busco_short_summary.txt')}
     """
 }
