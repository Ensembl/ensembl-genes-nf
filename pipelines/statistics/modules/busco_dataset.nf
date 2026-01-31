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
    val meta

    output:
    tuple val(meta), stdout, emit: busco_dataset_output
    path "versions.yml", emit: versions_file

    script:
    """
    echo "DEBUG: meta.core=${meta.dbname}, meta.species_id=${meta.species_id}, meta.taxon_id=${meta.taxon_id}" >&2

    if [[ !"${meta.busco_dataset}" ]]; then
    clade_selector.py -d ${params.busco_datasets_file} -t ${meta.taxon_id}  
    
    else 
    echo "${meta.busco_dataset}"
    fi
    # Create versions file
    PYTHON_VERSION=\$(python --version 2>&1 | sed 's/Python //')
    CLADE_VERSION=\$(clade_selector.py --version 2>&1 || echo "unknown")
    
    echo '"BUSCO_DATASET":' > versions.yml
    echo "  python: \$PYTHON_VERSION" >> versions.yml
    echo "  clade_selector: \$CLADE_VERSION" >> versions.yml
    """
}



