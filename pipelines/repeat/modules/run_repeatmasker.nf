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
process RUN_REPEATMASKER {
    label "python"
    tag "${meta.gca}:genome"

    // Multiple publishDir statements as needed
    //publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa", mode: "move"
    //publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.cat", mode: "move"
    //publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.masked", mode: "move"
    //publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.ori.out", mode: "move"
    //publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.out", mode: "move"
    //publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.tbl", mode: "move"
    //publishDir "${params.outDir}/repeatmasker/", pattern: "*.fa.rm.gtf", mode: "move"
    publishDir "${params.outDir}/repeatmasker/", pattern: "*.gtf", mode: "move"

    input:
    val meta

    output:
    val(meta), emit: repeatmasker_out
    path "versions.yml", emit: versions_file

    //path ("*.fa", emit: path_fasta),
    //path ("*.fa.cat", emit: path_fasta_cat),
    //path ("*.fa.masked", emit: path_fasta_masked),
    //path ("*.fa.ori.out", emit: path_fasta_ori_out),
    //path ("*.fa.out", emit: path_fa_out),
    //path ("*.fa.tbl", emit: path_fa_tbl),
    //path ("*.fa.rm.gtf", emit: path_fa_rm_gtf),
    //path ("*.gtf", emit: path_gtf)

    script:
    """
    run_repeatmasker --genome_file ${meta.genome_file} \
                    --output_dir ${params.outDir}/repeatmasker \
                    --repeatmasker_bin ${params.repeatmasker_path} \
                    --library ${params.outDir}/${meta.gca}/rm_library/${meta.gca}.repeatmodeler.fa \
                    --repeatmasker_engine ${params.engine_repeatmasker}
                    --num_threads ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -version 2>&1 | head -n 1 | sed 's/.*RepeatMasker version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

}

