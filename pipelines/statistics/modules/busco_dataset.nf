#!/usr/bin/env nextflow
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

process BUSCO_DATASET {

    label 'python'
    tag "${meta.core}"

    input:
    //tuple val(gca),  val(dbname), val(species_id)
    val meta

    output:
    //tuple val(gca), val(dbname), val(species_id), stdout
    tuple val(meta), stdout, emit: busco_dataset_output
    path "versions.yml", emit: versions_file
    

    
    script:
    """
    if [[ "${meta.taxon_id}" == "UNKNOWN" ]]; then
    TAXON_ID=\$(python utils.py \
    --db ${meta.core} \
    --key species.taxonomy_id \
    --species-id ${meta.species_id} \
    --host ${params.host} \
    --port ${params.port} \
    --user ${params.user_r}
    )
    else
    TAXON_ID="${meta.taxon_id}"
    fi

    clade_selector.py -d ${params.busco_datasets_file} -t "\$TAXON_ID"
    
    cat <<EOF > versions.yml
    BUSCO_DATASET:
      python: \$(python --version 2>&1 | sed 's/Python //')
        clade_selector: \$(clade_selector.py --version 2>&1 || echo "unknown")
    EOF
    """
}



