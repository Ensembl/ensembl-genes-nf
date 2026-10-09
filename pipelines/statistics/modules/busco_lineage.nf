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
/*BUSCO_LINEAGE process to run BUSCO in genome or protein mode
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca, busco_mode
- fasta_file: path to the genome or protein FASTA file
Outputs:
- busco_lineage_output: tuple of meta and BUSCO summary files
- versions_busco_<mode>.yml: versions file containing software versions used
*/

process BUSCO_LINEAGE {
    label 'busco'
    tag "${meta.gca}:busco_${meta.busco_mode}"

    publishDir { "${params.outdir}/${meta.gca}" }, mode: 'copy',
        saveAs: { filename -> filename.startsWith("busco_${meta.busco_mode}") ? filename : null }
    publishDir { "${params.outdir}/${meta.gca}" }, mode: 'copy', pattern: { "versions_busco_${meta.busco_mode}.yml" }
    afterScript "sleep ${params.files_latency}"
    maxForks 10

    input:
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path("busco_${meta.busco_mode}/*.txt"), emit: busco_lineage_output
    path "versions_busco_${meta.busco_mode}.yml", emit: versions_file
    path "busco_${meta.busco_mode}", emit: busco_full_output

    script:
    def busco_mode_arg = meta.busco_mode == 'protein' ? 'proteins' : 'genome'
    def busco_outdir = "busco_${meta.busco_mode}"
    def versions_file = "versions_busco_${meta.busco_mode}.yml"
    log.info("Selected BUSCO dataset: ${meta.busco_dataset}")

    """
    busco -f \
        -i ${fasta_file} \
        --mode ${busco_mode_arg} \
        -l ${meta.busco_dataset} \
        -c ${task.cpus} \
        --out ${busco_outdir} \
        --offline \
        --download_path ${params.download_path}

    cat <<EOF > ${versions_file}
    "${task.process}:${meta.busco_mode}":
        busco: \$(busco --version 2>&1 | head -n1 | awk '{print \$2}')
    EOF
    """
}
