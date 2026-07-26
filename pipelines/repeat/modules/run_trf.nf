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
This process runs the TRF (Tandem Repeats Finder) tool to identify 
tandem repeats in a genome file.
It uses the run_trf script to perform the analysis and generates a 
GTF file containing the identified tandem repeats.
The output GTF file is saved in the "trf" directory under the output 
directory for the given GCA accession.
The process also generates a versions.yml file containing the version 
of TRF used.
*/
process RUN_TRF {
    label "python"
    tag "${meta.gca}:genome"

    publishDir "${params.outdir}/${meta.gca}/trf/", pattern: "**/*.gtf", mode: "copy"

    input:
    val(meta)

    output:
    tuple val(meta), path("*.gtf"), emit: trf_out
    path "versions.yml", emit: versions_file
    
    script:
    """
    run_trf --genome_file ${meta.genome_file} \
                    --output_dir . \
                    --trf_bin ${params.trf_path} \
                    --match_score ${params.trf_match_score} \
                    --mismatch_score ${params.trf_mismatch_score} \
                    --delta ${params.trf_delta} \
                    --pm ${params.trf_pm} \
                    --pi ${params.trf_pi} \
                    --minscore ${params.trf_minscore} \
                    --maxperiod ${params.trf_maxperiod} \
                    --num_threads ${task.cpus} \
                    --bedtools_bin /opt/linuxbrew/bin/bedtools
    mv trf_output/annotation.gtf ${meta.gca}_trf.gtf                
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trf: \$(trf -version 2>&1 | head -n 1 | sed 's/.*Tandem Repeats Finder version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

}
