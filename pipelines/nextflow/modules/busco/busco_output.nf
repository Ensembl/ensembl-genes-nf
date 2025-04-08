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

process BUSCO_OUTPUT {
    label 'local'
    tag "${organism_name}:${insdc_acc}"
    publishDir "${params.outDir}/${publish_dir_name}/statistics", mode: 'copy'
    
    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source),
            path(summary_file)
        val(datatype)

    output:
        tuple val(insdc_acc), val(dbname), val(formated_sci_name), 
            val(publish_dir_name), path("*_short_summary.txt")
    
    script:
        formated_sci_name = organism_name.replaceAll("\\s", "_")
        publish_dir_name = formated_sci_name + '/' + insdc_acc + '/' + annotation_source
        species_lc = formated_sci_name.toLowerCase()

        if (dbname =='dummycore') {
            publish_dir_name = insdc_acc
            species_lc = 'species=NA'
        }

        if (datatype == "genome") {
            busco_sum = "genome_busco"
        }
        else if (datatype == "protein") {
            busco_sum = "protein_busco"
        }

        gca_string = insdc_acc.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")
        summary_name = [species_lc, gca_string, busco_sum, "short_summary.txt"].join("_")
        """
        sed '/Summarized benchmarking in BUSCO notation for file/d' $summary_file > $summary_name
        """
}
