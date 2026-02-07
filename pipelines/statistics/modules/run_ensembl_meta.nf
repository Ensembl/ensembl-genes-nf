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
/*RUN_ENSEMBL_META process to fetch Ensembl core database metadata SQL files
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca   
Outputs:
- ensembl_meta_output: tuple of meta and generated SQL files
- versions.yml: versions file containing software versions used
*/
process RUN_ENSEMBL_META {
    label 'python'
    tag "${meta.gca}"
    publishDir "${params.outdir}/${meta.gca}", mode: 'copy'
    afterScript "sleep ${params.files_latency}"

    input:
    val meta

    output:
    tuple val(meta), path("*.sql"), emit: ensembl_meta_output
    path "versions.yml", emit: versions_file

    script:
    """
    python ${params.enscode}/ensembl-genes/src/python/ensembl/genes/metadata/core_meta_data.py \
    --output_dir core_statistics --db_name ${meta.dbname} \
    --host ${params.host} --port ${params.port}  \
    --team ${params.team}  \
    --production_name ${meta.production_name}
    ln -s core_statistics/*.sql .
    # Create versions file
    PYTHON_VERSION=\$(python --version 2>&1 | awk '{print \$2}')
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$PYTHON_VERSION
    END_VERSIONS
    """
}
