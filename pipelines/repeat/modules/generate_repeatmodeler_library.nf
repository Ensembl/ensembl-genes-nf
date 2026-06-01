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
    tag "$meta.gca:run_repeatmodeler"
    label 'repeatmodeler'
    publishDir "${params.outdir}/${meta.gca}/library", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    tuple val(meta), path(genome_file)

    output:
    tuple val(meta), path("*-families.fa"), path("*-families.stk"), path("*-rmod.log"), emit: repeatmodeler_library_out
    path "versions.yml", emit: versions_file
    script:
    """
    echo "Running RepeatModeler for ${meta.gca} using genome file ${genome_file}"
    ${params.builddatabase_path} -name ${meta.gca}.repeatmodeler  ${genome_file}
    RepeatModeler -engine ${params.engine_repeatmodeler} -threads ${task.cpus} -database ${meta.gca}.repeatmodeler
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version 2>&1 | sed -n 's/.*RepeatModeler version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
