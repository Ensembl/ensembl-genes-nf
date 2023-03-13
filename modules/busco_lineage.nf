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


process BUSCO_LINEAGE {
  //label 'ncbi_taxonomy'
  label 'default'
  input:
  tuple val(gca), val(assembly_name), val(species_name)
  storeDir "${params.outDir}/${gca}/data/"
  //containerOptions "-B ${workDir}/:/busco_wd"
  output:
  val(gca), emit:gca
  val(assembly_name), emit: assembly_name
  //stdout emit: busco_lineage
  path "*.txt", emit: busco_lineage

  script:
  """
  mkdir -p ${params.outDir}/${gca}
  singularity run -H ${params.outDir}/${gca}:/home --bind ${params.outDir}/${gca}/:/data/:rw "${params.ncbi_taxonomy_singularity_path}" python /ncbi_taxa/get_busco_dataset.py --species_name "${species_name.trim()}"
  
  """
}
