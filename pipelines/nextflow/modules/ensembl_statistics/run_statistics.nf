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

process RUN_STATISTICS {
    label 'fetch_file'
    tag "core_statistics:${insdc_acc}"
    publishDir "${params.outDir}/${publish_dir_name}/", mode: 'copy'
    afterScript "sleep ${params.files_latency}"  // Needed because of file system latency
    maxForks 20    

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source)

    output:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), path("core_statistics/*.sql")

    script:
        def formated_sci_name = organism_name.replaceAll("\\s", "_")
        publish_dir_name = formated_sci_name + '/' + insdc_acc + '/' + annotation_source
        def corestats_outdir = "core_statistics"
        def stats_script = file("${params.enscode}/ensembl-genes/src/perl/ensembl/genes/generate_species_homepage_stats.pl")
        """
        perl ${stats_script} \
            -dbname ${dbname} \
            -host ${params.host} \
            -port ${params.port} \
            -production_name ${production_name} \
            -output_dir ${corestats_outdir}
        """
}
