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

process BUSCO_LINEAGES {
    label 'busco'
    tag "${organism_name}:${insdc_acc}"
    storeDir "${params.cacheDir}/${insdc_acc}/busco_${busco_mode}/"
    afterScript "sleep ${params.files_latency}"  // Needed because of file system latency
    maxForks 10

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), 
            val(ortho_db), path (aa_or_genome_seqs)
        val (input_busco_mode)

    output:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source),
            path("${outdir}/*.txt"), emit: busco_report_output

    script:
        if ( input_busco_mode == 'protein' ) {
            def outdir = 'protein_output'
            def busco_mode = 'proteins'
        }
        else if( input_busco_mode == 'genome' ) {
            def outdir = 'genome_output'
            def busco_mode = input_busco_mode
        }
        else{
            error "Invalid alignment mode: ${input_busco_mode}"
        }

        """
        mkdir ${outdir}
        busco -f \
            -i ${aa_or_genome_seqs} \
            --mode ${busco_mode} \
            -l ${ortho_db} \
            -c ${task.cpus} \
            --out ${outdir} \
            --offline \
            --download_path ${params.download_path}
        """
}