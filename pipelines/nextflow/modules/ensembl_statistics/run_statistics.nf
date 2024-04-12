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
    label 'omamer'
    tag "$db.species:$db.gca"
    storeDir "$cache_dir/$gca/omark_output"
    afterScript "sleep $params.files_latency"  // Needed because of file system latency
    input:
    file omamer_file

    output:
    path("proteins_detailed_summary.txt"), emit: summary_file

    script:
    """
    omark -f ${omamer_file} -d ${params.omamer_database} -o omark_output
    """
}