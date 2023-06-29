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
    // rename busco summary file in <production name>_gca_busco_short_summary.txt
    label 'default'
    publishDir "$output_dir/${db.species}/${db.gca}/statistics",  mode: 'copy'

    input:
    tuple val(db), path(summary_file, stageAs: "short_summary.txt")
    val(output_dir)
    val(name)

    output:
    path("*busco_short_summary.txt")

    script:
    def species = db.species.toLowerCase()
    def gca = db.gca.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")
    def summary_name = [species, gca, name, "short_summary.txt"].join("_")
    """
    sed '/Summarized benchmarking in BUSCO notation for file/d' $summary_file > $summary_name
    """
}
