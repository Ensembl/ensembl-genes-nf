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
/*BUSCO_GENOME_LINEAGE process to run BUSCO in genome mode
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca
- genome_file: path to the genome FASTA file
Outputs:
- busco_genome_lineage_output: tuple of meta and BUSCO output files
- versions_busco_genome.yml: versions file containing software versions used
*/

process BUSCO_GENOME_LINEAGE {
    label "busco"
    tag "${meta.gca}:busco_genome"
    //publishDir "${params.outdir}/${meta.gca}", mode: 'copy', pattern: "busco_genome/*"
    publishDir "${params.outdir}/${meta.gca}", mode: 'copy',
    saveAs: { filename -> filename.startsWith('busco_genome') ? filename : null }
    publishDir "${params.cacheDir}/${meta.gca}/busco_genome", mode: 'copy', pattern: "versions_busco_genome.yml"
    afterScript "sleep ${params.files_latency}"
    // Needed because of file system latency


    input:
    tuple val(meta), path(genome_file)

    output:
    tuple val(meta), path("busco_genome/*.txt"), emit: busco_genome_lineage_output
    path "versions_busco_genome.yml", emit: versions_file
    path "busco_genome", emit: busco_full_output

    script:
    log.info("Selected BUSCO dataset: ${meta.busco_dataset}")

    """
    busco -f \
    -i ${genome_file} \
    --mode genome \
    -l ${meta.busco_dataset} \
    -c ${task.cpus} \
    --out busco_genome \
    --offline \
    --download_path ${params.download_path}


    cat <<EOF > versions_busco_genome.yml
    BUSCO_GENOME_LINEAGE:
        busco: \$(busco --version 2>&1 | head -n1 | awk '{print \$2}')
    EOF
    """
}
