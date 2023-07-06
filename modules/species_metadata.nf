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

// Get species name, gca accession and annotation source from meta table
process SPECIES_METADATA {
  label 'default'
  tag "$dbname"

  input:
  val dbname
  val output_dir
  val project

  output:
  stdout

  script:
  """
  function get_metadata {
    KEY=\$1
    mysql -N -u ${params.user} \
             -h ${params.host} \
             -P ${params.port} \
             -D $dbname \
             -e "SELECT meta_value FROM meta WHERE meta_key='\$KEY'"
  }
  SPECIES=\$(get_metadata "species.scientific_name" | sed 's/ /_/g')
  GCA=\$(get_metadata "assembly.accession")
  SOURCE=\$(get_metadata "species.annotation_source")
  if [ "\$SOURCE" = "" ]; then SOURCE=""; fi
  BRC_COMPONENT=\$(get_metadata "BRC4.component")
  if [ "\$BRC_COMPONENT" = "" ]; then BRC_COMPONENT=""; fi
  BRC_ORGANISM=\$(get_metadata "BRC4.organism_abbrev")
  if [ "\$BRC_ORGANISM" = "" ]; then BRC_ORGANISM=""; fi

  PUBLISH_DIR=$output_dir
  if [ $project == 'ensembl' ]; then
    PUBLISH_DIR="$output_dir/\$SPECIES/\$GCA"
  fi
  if [ $project == 'BRC' ]; then
    PUBLISH_DIR="$output_dir/\$BRC_COMPONENT/\$BRC_ORGANISM"
  fi

  echo "name,species,gca,publish_dir,source,brc_component,brc_organism"
  echo "$dbname,\$SPECIES,\$GCA,\$PUBLISH_DIR,\$SOURCE,\$BRC_COMPONENT,\$BRC_ORGANISM"
  """
}
