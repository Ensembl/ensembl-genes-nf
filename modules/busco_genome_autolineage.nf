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

// run Busco in genome mode 
process BUSCO_GENOME_AUTOLINEAGE {

  label 'busco'

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
}