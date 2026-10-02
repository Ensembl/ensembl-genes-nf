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
/*BUSCO_ASSEMBLY_DB process to load BUSCO genome results into the assembly metadata database
Inputs:
- meta: metadata map containing gca, busco_mode
- summary_file: path to the BUSCO genome summary file
- asm_metadata: JSON string with the assembly metadata DB connection params
Outputs:
- *_busco_genome_asmdb.json: JSON file containing the BUSCO metrics loaded
- versions_busco_asmdb.yml: versions file containing software versions used
*/
process BUSCO_ASSEMBLY_DB {

    label 'python'
    tag "${meta.gca}"
    cache false
    publishDir "${params.outdir}/${meta.gca}", mode: 'copy'
    afterScript "sleep ${params.files_latency}"

    input:
    tuple val(meta), path(summary_file)
    val asm_metadata

    output:
    tuple val(meta), path("${meta.gca}_busco_genome_asmdb.json"), emit: metrics_json
    path "versions_busco_asmdb.yml", emit: versions_file

    when:
    params.update_registry

    script:
    // Single-quote the JSON for the shell, escaping any single quotes it contains
    def asm_metadata_arg = "'" + asm_metadata.replace("'", "'\\''") + "'"
    """
    busco_metakeys_to_asmdb.py \
        --gca ${meta.gca} \
        --file ${summary_file} \
        --asm-metadata ${asm_metadata_arg} \
        --output-json ${meta.gca}_busco_genome_asmdb.json

    cat <<-END_VERSIONS > versions_busco_asmdb.yml
    "${task.process}":
        busco_metakeys_to_asmdb.py: \$(busco_metakeys_to_asmdb.py --version 2>&1 | grep -oP 'version \\K[0-9.]+' || echo "unknown")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
