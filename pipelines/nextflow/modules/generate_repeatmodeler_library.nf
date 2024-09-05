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

process GENERATE_REPEATMODELER_LIBRARY {
    tag "$gca:run_repeatmodeler"
    label 'repeatmodeler'
    publishDir "${params.outDir}/${gca}/rm_database", mode: 'copy'
    publishDir "${params.outDir}/${gca}/repeatmodeler_output", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    maxForks 10

    input:
    tuple val(gca)

    output:
    tuple val(gca), path("${gca}*.repeatmodeler_db"), path ("repeatmodeler_output/*")

    script:
    """
    ${params. builddatabase_path} -name ${gca}.repeatmodeler -engine  ${params.engine_repeatmodeler} ${params.outDir}/${gca}/ncbi_dataset/${gca}*.fna
    ${params.repeatmodeler_path} -engine ${params.engine_repeatmodeler} -pa 10 -database ${params.outDir}/${gca}/rm_database/${gca}.repeatmodeler -dir ${params.outDir}/${gca}/repeatmodeler_output
    mv ${params.outDir}/${gca}/repeatmodeler_output/${gca}.repeatmodeler-families.fa ${params.outDir}/${gca}/repeatmodeler_output/${gca}.repeatmodeler.fa
    cp ${params.outDir}/${gca}/repeatmodeler_output/${gca}.repeatmodeler.fa ${params.outDir}/${gca}/rm_library
    """
}
