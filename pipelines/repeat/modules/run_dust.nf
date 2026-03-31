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
process RUN_DUST {
    label "python"
    tag "${meta.gca}:genome"


    publishDir "${params.outDir}/dust/", pattern: "*.gtf", mode: "move"

    input:
    val meta

    output:
    val(meta), emit: dust_out
    path "versions.yml", emit: versions_file

    script:
    """
    run_dust --genome_file ${meta.genome_file} \
                    --output_dir ${params.outDir}/dust \
                    --dust_bin ${params.dust_path} \
                    --num_threads ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dust: \$(Dust -version 2>&1 | head -n 1 | sed 's/.*Dust version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

}

