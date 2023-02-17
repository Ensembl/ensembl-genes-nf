def get_species_name(name) {
  final m = name.tr('[A-Z]','[a-z]').tr('.','v').split('/')[0]
  return m
}

def get_gca(name) {
  final m = name.tr('[A-Z]', '[a-z]').tr('.', 'v').replaceAll("_","").split('/').getAt(1)
  return m
}

def concatString(string1, string2, string3){
 return string1 + '_'+string2 + '_'+string3
}

/* Get species name and accession from meta table to build the output directory tree */
process SPECIESOUTDIR {
  cpus 1
  memory { 2.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  val db
  output:
  tuple stdout, val(db)
  script:
  // get <Production name>/GCA
  """
  mysql -N -u ${params.user}  -h ${params.host} -P ${params.port} -D $db < "${params.meta_query_file}"
  """
}

/* Dump canonical translations */
process FETCHPROTEINS {
  cpus 1
  memory { 9.GB * task.attempt }
  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  tuple val(species_dir),val(db)
  
  storeDir "${params.outDir}/${species_dir.trim()}/fasta/"

  output:
  path "translations.fa", emit: fasta
  val "${species_dir.trim()}", emit: output_dir
  val db, emit:db_name

  beforeScript "export ENSCODE=${params.enscode}"
  beforeScript "source $ENSCODE/ensembl-genes-nf/supplementary_files/perl5lib.sh"
  //beforeScript "ENSCODE=${params.enscode} source ${projectDir}/supplementary_files/perl5lib.sh"
  script:
  """
  perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl -host ${params.host} -port ${params.port} -dbname $db -user ${params.user} -dnadbhost ${params.host} -dnadbport ${params.port} -dnadbname $db -dnadbuser ${params.user} -canonical_only 1 -file translations.fa  ${params.dump_params}
  """
}


process CREATEOMAMER {

  cpus 20
  memory { 10.GB * task.attempt }

  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  file translations
  val outdir
  val db

  output:
  path "*.omamer", emit: omamer_file
  val outdir, emit:species_outdir

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.outDir}/${outdir}/",  mode: 'copy'

  script:
  """
  singularity exec /hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/omark.sif omamer search --query ${translations} --db ${params.omamer_database} --score sensitive --out proteins.omamer
  """
}


process RUNOMARK {

  cpus 20
  memory { 20.GB * task.attempt }

  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
  file omamer_file
  val outdir
  

  output:
  path("omark_output/*"), emit: summary_file
  val outdir, emit:species_outdir

  publishDir "${params.outDir}/${outdir}/",  mode: 'copy'

  script:
  """
  singularity exec /hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/omark.sif omark -f ${omamer_file} -d ${params.omamer_database} -o omark_output
  """
}
