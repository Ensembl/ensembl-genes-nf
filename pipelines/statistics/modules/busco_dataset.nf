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
/*BUSCO_DATASET process to select BUSCO dataset based on taxon_id
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca
Outputs:
- busco_dataset_output: tuple of meta and selected busco_dataset
- versions.yml: versions file containing software versions used
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
    def busco_dataset = params.busco_dataset ?: meta.busco_dataset ?: ''
    busco_dataset = busco_dataset ? busco_dataset.trim() : ''
    """
    export PYTHONPATH="${params.enscode}/ensembl-genes/src/python:\${PYTHONPATH:-}"

    if [[ -z "${busco_dataset}" ]]; then
        python ${params.enscode}/ensembl-genes/src/python/ensembl/genes/metrics/busco_lineage_selector.py -d ${params.busco_datasets_file} -t ${meta.taxon_id}
    else
        echo "${busco_dataset}"
    fi
    # Create versions file
    PYTHON_VERSION=\$(python --version 2>&1 | sed 's/Python //')

    
    echo '"BUSCO_DATASET":' > versions.yml
    echo "  python: \$PYTHON_VERSION" >> versions.yml
    """
}
