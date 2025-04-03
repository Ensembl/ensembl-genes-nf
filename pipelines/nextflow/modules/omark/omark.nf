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

process OMARK {
    label 'omamer'
    tag "${organism_name}:${insdc_acc}"
    publishDir "${params.outDir}/${publish_dir_name}/", mode: 'copy'

    container "${params.omark_singularity_path}"

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), path(omamer_file)

    output:
        tuple val(insdc_acc), val(dbname), val(formated_sci_name), val(publish_dir_name), path("omark_output/*_summary.txt"), path("omark_output/*")

    script:
        formated_sci_name = organism_name.replaceAll("\\s", "_")
        publish_dir_name = formated_sci_name + '/' + insdc_acc + '/' + annotation_source
    
        """
        omark -f ${omamer_file} -d ${params.omamer_database} -o omark_output
        """
}
