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
    tag "$db.gca:protein"
    label 'fetch_file'
    storeDir "$cache_dir/${db.gca}/fasta/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    val(db)
    val cache_dir

    output:
    //tuple val(db), val(busco_dataset), path("*_translations.fa")
    path("*_translations.fa"), emit:fasta
    script:
    def translations_file = "${db.name}_translations.fa"
    """
    perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl \
        -host ${params.host} \
        -port ${params.port} \
        -dbname ${db.name} \
        -user ${params.user} \
        -file $translations_file \
        ${params.dump_params}
    """
}
