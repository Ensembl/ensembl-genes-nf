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
This process runs RepeatModeler on a given genome file to generate a RepeatModeler library. 
It uses the BuildDatabase tool to create a database from the genome file and then runs 
RepeatModeler with the specified engine and number of threads. 
The output files are renamed to include the GCA accession in their names 
and are saved in the "library" directory under the output directory for 
the given GCA accession. The process also generates a versions.yml file 
containing the version of RepeatModeler used.
*/

process GENERATE_REPEATMODELER_LIBRARY {
    tag "$meta.gca:run_repeatmodeler"
    label 'repeatmodeler'
    publishDir "${params.outdir}/${meta.gca}/library", mode: 'copy'
    afterScript "sleep $params.files_latency"  // Needed because of file system latency

    input:
    val(meta)

    output:
    tuple val(meta), path("*repeatmodeler.fa"), path("*families.stk.gz"), path("*rmod.log"), emit: repeatmodeler_library_out
    path "versions.yml", emit: versions_file
    script:
    """
    echo "Running RepeatModeler for ${meta.gca} using genome file ${meta.genome_file}"
    BuildDatabase -name ${meta.gca}.repeatmodeler  ${meta.genome_file}
    echo "Database files after BuildDatabase:"
    RepeatModeler -engine ${params.engine_repeatmodeler} -threads ${task.cpus} -database ${meta.gca}.repeatmodeler
    # Rename outputs to the desired published names
    mv ${meta.gca}.repeatmodeler-families.fa ${meta.gca}.repeatmodeler.fa
    gzip -f ${meta.gca}.repeatmodeler-families.stk
    mv ${meta.gca}.repeatmodeler-families.stk.gz ${meta.gca}.families.stk.gz
    mv ${meta.gca}.repeatmodeler-rmod.log ${meta.gca}.rmod.log
                    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version 2>&1 | sed -n 's/.*RepeatModeler version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
