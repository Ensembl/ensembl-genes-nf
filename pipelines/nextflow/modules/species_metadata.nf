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

  function has_any {
    TABLE=\$1
    mysql -N -u ${params.user} \
            -h ${params.host} \
            -P ${params.port} \
            -D $dbname \
            -e "SELECT stable_id FROM \$TABLE LIMIT 1" | wc -l
  }

  SPECIES=\$(get_metadata "species.scientific_name" | sed 's/ /_/g')
  GCA=\$(get_metadata "assembly.accession")
  TAXON_ID=\$(get_metadata "species.taxonomy_id")
  PRODUCTION_NAME=\$(get_metadata "species.production_name")
  SOURCE=\$(get_metadata "species.annotation_source")
  if [ "\$SOURCE" = "" ]; then SOURCE=""; fi
  BRC_COMPONENT=\$(get_metadata "BRC4.component")
  if [ "\$BRC_COMPONENT" = "" ]; then BRC_COMPONENT=""; fi
  BRC_ORGANISM=\$(get_metadata "BRC4.organism_abbrev")
  if [ "\$BRC_ORGANISM" = "" ]; then BRC_ORGANISM=""; fi
  HAS_GENES=\$(has_any "gene")

  PUBLISH_DIR=${params.outDir}
  if [ ${params.project} == 'ensembl' ]; then
    PUBLISH_DIR="${params.outDir}/\$SPECIES/\$GCA/\$SOURCE/"
  fi
  if [ ${params.project} == 'brc' ]; then
    PUBLISH_DIR="$output_dir/\$BRC_COMPONENT/\$BRC_ORGANISM"
  fi

  echo "name,species,gca,taxon_id,production_name,publish_dir,source,has_genes,brc_component,brc_organism"
  echo "$dbname,\$SPECIES,\$GCA,\$TAXON_ID,\$PRODUCTION_NAME,\$PUBLISH_DIR,\$SOURCE,\$HAS_GENES,\$BRC_COMPONENT,\$BRC_ORGANISM"
  """
}