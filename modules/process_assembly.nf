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

// Import utility functions
include { build_ncbi_path } from './utils.nf'
include { concatString } from './utils.nf'

process PROCESS_ASSEMBLY {
  //memory { 8.GB * task.attempt }
  label 'default'

  input:
  //tuple val(gca), val(assembly_name)
  val gca
  val assembly_name
  val busco_lineage
  storeDir "${params.outDir}/${gca}/data/"
  
  output:
  val(gca), emit:gca
  path "*.fna", emit: genome_file
  val busco_lineage, emit:busco_lineage

  script:
  """
  wget  ${build_ncbi_path("${gca}", "${assembly_name}")}
  gzip -d -f ${concatString("${gca}", "${assembly_name.replaceAll(" ","_")}", 'genomic.fna.gz')}
  """
}
