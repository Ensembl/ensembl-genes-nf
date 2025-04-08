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
    tag "${organism_name}:${insdc_acc}"
    label 'fetch_file'
    storeDir "${params.cacheDir}/${insdc_acc}/fasta/"
    afterScript "sleep ${params.files_latency}"  // Needed because of file system latency
    maxForks 20

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), val(ortho_db)

    output:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source),
            val(ortho_db), path("*_translations.fa"), emit: translations

    script:
        def translations_file = production_name +"_translations.fa"
        def dump_translations_script = file("${params.enscode}/ensembl-analysis/scripts/protein/dump_translations.pl")
        """
        perl ${dump_translations_script} \
            -host ${params.host} \
            -port ${params.port} \
            -dbname ${dbname} \
            -user ${params.user_r} \
            -file ${translations_file} \
            ${params.canonical_only}
        """
}
