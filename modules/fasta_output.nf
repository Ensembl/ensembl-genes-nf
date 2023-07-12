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


include { make_publish_dir } from './utils.nf'


// dump unmasked dna sequences from core db 
process FASTA_OUTPUT {
  tag "$db.species"
  label "default"
  publishDir { make_publish_dir(db.publish_dir, project, name) },  mode: 'copy'

  input:
  tuple val(db), val(busco_dataset), path(fasta_file)
  val project
  val name

  output:
  path(fasta_file)

  script:
  """
  echo 'Using publish_dir'
  """
}
