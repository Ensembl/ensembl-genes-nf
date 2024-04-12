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

//include { make_publish_dir } from '../utils.nf'

process BUSCO_OUTPUT {
    label 'default'
    tag "busco_output"
    publishDir "$publish_dir/statistics", mode: 'copy'

    input:
    tuple val(db), path(summary_file), val(datatype)
    def publish_dir = db ? db.publish_dir : $params.outDir+"/"+db.gca
    

    output:
    path("*_short_summary.txt"), emit:summary_file

    script:
    def summary_name = summary_file
    def name = ""
    if (datatype == "genome") {
        name = "genome_busco"
    } else if (datatype == "protein") {
            name = "protein_busco"
    }
    if ($params.project == 'brc') {
        summary_name = [name, "short_summary.txt"].join("_")
    }
    else {
        def species = db ? db.species?.toLowerCase() : "species"
        def gca = db.gca.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")
        summary_name = [species, gca, name, "short_summary.txt"].join("_")
    } 
    """
    sed '/Summarized benchmarking in BUSCO notation for file/d' $summary_file > $summary_name
    """
}
