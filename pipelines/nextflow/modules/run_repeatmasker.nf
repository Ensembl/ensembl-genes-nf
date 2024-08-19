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
    label 'run_repeatmasker'
    tag "$gca:genome"
    publishDir "${params.outDir}/repeatmasker/", pattern:"*.fa" , mode "copy"
    publishDir "${params.outDir}/repeatmasker/", pattern:"*.fa.cat" , mode "copy"
    publishDir "${params.outDir}/repeatmasker/", pattern:"*.fa.masked" , mode "copy"
    publishDir "${params.outDir}/repeatmasker/", pattern:"*.fa.ori.out" , mode "copy"
    publishDir "${params.outDir}/repeatmasker/", pattern:"*.fa.out" , mode "copy"
    publishDir "${params.outDir}/repeatmasker/", pattern:"*.fa.tbl" , mode "copy"
    publishDir "${params.outDir}/repeatmasker/", pattern:"*.fa.rm.gtf" , mode "copy"

    maxForks 20

    input:
     tuple val(gca)

    output:
     tuple val(gca),
     path "*.fa", emit: slice_fasta
     path "*.fa.cat", emit slice_fasta_cat
     path "*.fa.masked", emit slice_fasta_masked
     path "*.fa.ori.out", emit slice_fasta_ori_out
     path "*.fa.out", emit slice_fa_out
     path "*.fa.tbl", emit slice_fa_tbl
     path "*.fa.rm.gtf", emit slice_fa_rm_gtf

    script:
    """
    chmod +x $projectDir/bin/repeatmasker.py #set executable permissions

    repeatmasker.py --genome_file ${params.cacheDir}/${gca}/ncbi_dataset/${gca}*.fna --output_dir ${params.outDir}/repeatmasker --repeatmasker_bin ${params.repeatmasker_path} --library ${params.cacheDir}/${gca}/rm_library/${gca}.repeatmodeler.fa --repeatmasker_engine ${params.engine}
    """
}
