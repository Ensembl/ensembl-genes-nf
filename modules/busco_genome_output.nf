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

// Import utility functions
include { concatString } from './utils.nf'

process BUSCO_GENOME_OUTPUT {
    //rename busco summary file in <production name>_gca_genome_busco_short_summary.txt
    label 'default'
    publishDir "$output_dir/${db.species}/${db.gca}/genome",  mode: 'copy'

    input:
    tuple val(db), path(summary_file)
    val(output_dir)

    output:
    path("statistics")

    script:
    def stats_dir = "statistics"
    def summary_name = concatString(db.species, db.gca, 'genome_busco_short_summary.txt')
    """
    mkdir $stats_dir
    sed '/Summarized benchmarking in BUSCO notation for file/d' $summary_file > $stats_dir/$summary_name
    """
}
