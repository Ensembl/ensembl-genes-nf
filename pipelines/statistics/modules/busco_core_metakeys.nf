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
/*BUSCO_CORE_METAKEYS process to add BUSCO results to core database metadata
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca
- summary_file: path to the BUSCO summary file  
Outputs:
- versions.yml: versions file containing software versions used

*/
process BUSCO_CORE_METAKEYS {

    label 'python'
    tag "${meta.gca}"
    publishDir "${params.outdir}/${meta.gca}", mode: 'copy'
    afterScript "sleep ${params.files_latency}"

    input:
    tuple val(meta), path(summary_file)

    output:
    path "versions.yml", emit: versions_file, optional: true

    when:
    params.apply_busco_metakeys

    script:

    """
    busco_metakeys_patch.py \
    -db ${meta.dbname} -file ${summary_file} \
    -output_dir "${params.outdir}/${meta.gca}/"  -host ${params.host} \
    -port ${params.port} -user ${params.user}  \
    -password ${params.password} -species_id ${meta.species_id} \
    -run_query true
    
    # Create versions file
    PYTHON_VERSION=\$(python --version 2>&1 | awk '{print \$2}')

    echo '"BUSCO_CORE_METAKEYS":' > versions.yml
    echo "  python: \$PYTHON_VERSION" >> versions.yml
    """
}
