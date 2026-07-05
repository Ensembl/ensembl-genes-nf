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
/*This process runs DustMasker to identify low-complexity regions in a genome file. 
It uses the dustmasker tool to perform the analysis and generates a GTF file 
containing the identified low-complexity regions. The output GTF file is saved 
in the "dust" directory under the output directory for the given GCA accession. 
The process also generates a versions.yml file containing the version of Dust used.*/

process RUN_DUST {
    label "python"
    tag "${meta.gca}:genome"

    publishDir "${params.outdir}/${meta.gca}/dust/", pattern: "**/*.gtf", mode: "copy"
    input:
    val(meta)

    output:
    tuple val(meta), path("**/*.gtf"), emit: dust_out
    path "versions.yml", emit: versions_file

    script:
    """
    run_dust --genome_file ${meta.genome_file} \
                    --output_dir . \
                    --dust_bin /opt/linuxbrew/bin/dustmasker \
                    --num_threads ${task.cpus}  \
                    --bedtools_bin /opt/linuxbrew/bin/bedtools
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dust: \$(Dust -version 2>&1 | head -n 1 | sed 's/.*Dust version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
