#!/usr/bin/env nextflow
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
def get_species_name(name) {
  final m = name.tr('[A-Z]','[a-z]').tr('.','v').split('/')[0]
  return m
}

def get_gca(name) {
  final m = name.tr('[A-Z]', '[a-z]').tr('.', 'v').replaceAll("_","").split('/').getAt(1)
  return m
}



def get_gca_path(gca, assembly_name) {
  final gca_splitted = gca.replaceAll("_","").tokenize(".")[0].split("(?<=\\G.{3})").join('/')
  //return gca_splitted
  return  'https://ftp.ncbi.nlm.nih.gov/genomes/all'  + '/' + gca_splitted + '/' + "$gca" +'_' + assembly_name.replaceAll(" ","_") + '/' + "$gca" + '_' + assembly_name.replaceAll(" ","_") + '_genomic.fna.gz'
  //return gca + '_' + assembly_name.replaceAll(" ","_") + '_genomic.fna.gz'    ?<=\\G.{3}
}


def concatString(string1, string2, string3){
 return string1 + '_'+string2 + '_'+string3
}


process PROCESSASSEMBLY {
  memory { 8.GB * task.attempt }
  
  input:
  tuple val(gca), val(assembly_name)

  storeDir "${params.outDir}/${gca}/data/"
  //publishDir "${params.outDir}/${gca}/data",  mode: 'copy'
  output:
  val(gca), emit:gca
  path "*.fna", emit: genome_file
  //print ${get_gca_path(${gca}, ${assembly_name})
  script:
  //"""
  //x= get_gca_path("${gca}", "${assembly_name}")
  //print x
  """
  wget  ${get_gca_path("${gca}", "${assembly_name}")}
  gzip -d -f ${concatString("${gca}", "${assembly_name.replaceAll(" ","_")}", 'genomic.fna.gz')}
  """
  //gzip -d ${concatString("${gca}",  "${assembly_name}.replaceAll(' ','_')", 'genomic.fna')}   ${get_gca_path("${gca}", "${assembly_name}")};
  //"""
  
}

/* run Busco in genome mode */
process BUSCOGENOME {

  cpus 40
  memory { 40.GB * task.attempt }
  time "24h"

  errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
  maxRetries 2
  module 'singularity-3.7.0-gcc-9.3.0-dp5ffrp'
  container "ezlabgva/busco:${params.busco_version}"
  containerOptions "-B ${params.outDir}/:/busco_wd"

  input:
  val gca
  file genome_file

  output:
  path "busco_output/*.txt", emit: summary_file
  publishDir "${params.outDir}/${gca}/",  mode: 'copy'

  script:
  """
  busco -f -i ${genome_file} -o busco_output --mode genome --auto-lineage -c ${task.cpus} 
  """
  //--offline --download_path ${params.download_path}
  //"""
}

