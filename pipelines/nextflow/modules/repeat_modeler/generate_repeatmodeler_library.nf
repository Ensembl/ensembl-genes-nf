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
    publishDir "${params.outDir}/${gca}/", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(species),val(gca),path(genome_file)

    output:
    tuple val(gca), val(species), path("*")

    script:
    """
    echo "Running RepeatModeler for ${gca} using genome file ${genome_file}"
    ${params.builddatabase_path} -name ${gca}.repeatmodeler -dir ${params.outDir}/${gca}
    singularity run ${params.repeatmodeler_path} RepeatModeler -engine ${params.engine_repeatmodeler} -threads ${task.cpus} -database ${gca}.repeatmodeler
    """
}
