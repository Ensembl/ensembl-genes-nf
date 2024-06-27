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

include { getMetaValue } from './utils.nf'

process FETCH_PROTEINS {
    tag "$gca:protein"
    label 'fetch_file'
    storeDir "${params.cacheDir}/$gca/fasta/"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 20

    input:
    tuple val(gca), val(dbname), val(busco_dataset)

    output:
    tuple val(gca), val(dbname), path("*_translations.fa"),val(busco_dataset) 

    script:
    scientific_name = getMetaValue(dbname, "species.production_name")[0].meta_value.toString().toLowerCase()
    translations_file = scientific_name +"_translations.fa"
    """
    perl ${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl \
        -host ${params.host} \
        -port ${params.port} \
        -dbname ${dbname} \
        -user ${params.user_r} \
        -file $translations_file \
        ${params.dump_params}
    """
}
