/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

process FETCH_GENOME {
  tag "$gca:genome"
  label 'fetch_file'
  storeDir "${params.cacheDir}/$gca/ncbi_dataset/"
  afterScript "sleep $params.files_latency"  // Needed because of file system latency
  maxForks 10

  input:
    tuple val(gca), val(dbname), val(busco_dataset)

  output:
    tuple val(gca), val(dbname), path("*.fa"), val(busco_dataset)
  
  script:
    genome_fasta = "genome_toplevel.fa"
    sequence_dumping_script = file("${params.enscode}/ensembl-analysis/scripts/sequence_dump.pl")

    """
    perl ${sequence_dumping_script} \
      -dbhost ${params.host} \
      -dbport ${params.port} \
      -dbname ${dbname} \
      -dbuser ${params.user_r} \
      -coord_system_name toplevel \
      -toplevel \
      -onefile \
      -nonref \
      -filename $genome_fasta
    """
}
