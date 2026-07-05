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
/*This process runs RepeatMasker on a given genome file using a 
specified RepeatModeler library. It uses the RepeatMasker tool
to perform the analysis and generates GTF files containing the identified 
repetitive regions. The output GTF files are saved in the "repeatmasker" 
directory under the output directory for the given GCA accession. 
The process also generates a versions.yml file containing the version 
of RepeatMasker used.*/
process RUN_REPEATMASKER {
    label "python"
    tag "${meta.gca}:genome"

    publishDir "${params.outdir}/${meta.gca}/repeatmasker/", pattern: "repeatmasker_output/*.gtf", mode: "copy"

    input:
    tuple val(meta), val(library_file)

    output:
    tuple val(meta), path("repeatmasker_output/*.gtf"), emit: repeatmasker_out
    path "versions.yml", emit: versions_file

    script:
    LD_LIBRARY_PATH='/hps/software/users/ensembl/genebuild/shared/libnsl/lib:/hps/software/users/ensembl/genebuild/shared/libnsl/libtirpc'
    """
    export LD_LIBRARY_PATH=/hps/software/users/ensembl/genebuild/shared/libnsl/lib:/hps/software/users/ensembl/genebuild/shared/libnsl/libtirpc:$LD_LIBRARY_PATH
    ls /hps/software/users/ensembl/genebuild/shared/libnsl/lib/libnsl.so.2
    /opt/linuxbrew/bin/RepeatMasker -help
    run_repeatmasker  --genome_file ${meta.genome_file} \
                    --output_dir . \
                    --repeatmasker_bin /opt/linuxbrew/bin/RepeatMasker \
                    --library ${library_file} \
                    --repeatmasker_engine ${params.engine_repeatmasker} \
                    --num_threads ${task.cpus} \
                    --bedtools_bin /opt/linuxbrew/bin/bedtools
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -version 2>&1 | head -n 1 | sed 's/.*RepeatMasker version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """

}
