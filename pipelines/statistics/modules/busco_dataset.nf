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
    tag "${meta.gca}"

    input:
    //tuple val(gca),  val(dbname), val(species_id)
    val meta

    output:
    tuple val(meta), stdout, emit: busco_dataset_output
    path "versions.yml", emit: versions_file

    script:
    """
    # Debug output to stderr (won't interfere with stdout capture)
    echo "DEBUG: meta.core=${meta.dbname}, meta.species_id=${meta.species_id}, meta.taxon_id=${meta.taxon_id}" >&2
    
#    if [[ "${meta.taxon_id}" == "UNKNOWN" ]]; then
#        echo "DEBUG: Fetching taxon_id from database..." >&2
#          TAXON_ID=\$(utils.py \
#            --db ${meta.dbname} \
#            --key species.taxonomy_id \
#            --species-id ${meta.species_id} \
#            --host ${params.host} \
#            --port ${params.port} \
#            --user ${params.user_r})
#        echo "DEBUG: Found TAXON_ID=\$TAXON_ID" >&2
#    else
#        TAXON_ID="${meta.taxon_id}"
#        echo "DEBUG: Using provided TAXON_ID=\$TAXON_ID" >&2
#    fi

    # Run clade_selector and output to stdout (this gets captured!)
#    echo "DEBUG: Running clade_selector with TAXON_ID=\$TAXON_ID" >&2
#    ${projectDir}/bin/clade_selector.py -d ${params.busco_datasets_file} -t "\$TAXON_ID"
    if [[ !"${params.busco_dataset}" ]]; then
    clade_selector.py -d ${params.busco_datasets_file} -t ${meta.taxon_id}  
    
    else 
    echo "${params.busco_dataset}"
    fi
    # Create versions file
    PYTHON_VERSION=\$(python --version 2>&1 | sed 's/Python //')
    CLADE_VERSION=\$(${projectDir}/bin/clade_selector.py --version 2>&1 || echo "unknown")
    
    echo '"BUSCO_DATASET":' > versions.yml
    echo "  python: \$PYTHON_VERSION" >> versions.yml
    echo "  clade_selector: \$CLADE_VERSION" >> versions.yml
    """
}



