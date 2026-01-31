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

process BUSCO_PROTEIN_LINEAGE {
    label 'busco'
    tag "$meta.gca"
    publishDir "${params.outdir}/${meta.gca}", mode: 'copy'
    publishDir "${params.outdir}/${meta.gca}", mode: 'copy', pattern: "versions_busco_protein.yml"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 10

    input:
    tuple val(meta), path(translation_file)

    output:
    tuple val(meta), path("busco_protein/*.txt"), emit: busco_protein_lineage_output
    path "versions_busco_protein.yml", emit: versions_file

    script:
    log.info("Selected BUSCO dataset: $meta.busco_dataset")

    """
    busco -f \
        -i ${translation_file} \
        --mode proteins \
        -l ${meta.busco_dataset} \
        -c ${task.cpus} \
        --out busco_protein \
        --offline \
        --download_path ${params.download_path}

    cat <<EOF > versions_busco_protein.yml
    BUSCO_PROTEIN_LINEAGE:
        busco: \$(busco --version 2>&1 | head -n1 | awk '{print \$2}')
    EOF
    """
}
