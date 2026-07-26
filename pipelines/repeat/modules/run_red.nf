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
/*
This process runs the Red to identify repetitive regions in a genome file. 
It uses the Red tool to perform the analysis and generates a GTF file 
containing the identified repetitive regions. The output GTF file is 
saved in the "red" directory under the output directory for the given 
GCA accession. The process also generates a versions.yml 
file containing the version of Red used.
*/
process RUN_RED {
    label "python"
    tag "${meta.gca}:genome"
    publishDir "${params.outdir}/${meta.gca}/red/", pattern: "**/*.gtf", mode: "copy"

    input:
    val(meta)

    output:
    tuple val(meta), path("*.gtf"), emit: red_out
    path "versions.yml", emit: versions_file    

    script:
    """
    run_red --genome_file ${meta.genome_file} \
                    --output_dir . \
                    --red_bin ${params.red_path}
    mv red_output/annotation.gtf ${meta.gca}_red.gtf  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        red: \$(Red --version 2>&1 | head -n 1 | sed 's/.*Red version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

}
