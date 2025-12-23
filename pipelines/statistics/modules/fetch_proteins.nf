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


process FETCH_PROTEINS {
    tag "$meta.gca:protein"
    label 'fetch_file'
    storeDir "${params.cacheDir}/$meta.gca/fasta/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 20

    input:
    //tuple val(gca), val(dbname), val(species_id), val(busco_dataset)
    val(meta)

    output:
    //tuple val(gca), val(dbname), path("*_translations.fa"),val(busco_dataset) , val(species_id)
    tuple val(meta), path("*translations.fa"), emit: protein_file_output
    path "versions.yml", emit: versions_file

    script:
    translations_file = "translations.fa"
    """
    perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl \
        -host ${params.host} \
        -port ${params.port} \
        -dbname ${meta.dbname} \
        -user ${params.user_r} \
        -file $translations_file \
        --species_id ${meta.species_id} \
        ${params.dump_params}
    """
}
