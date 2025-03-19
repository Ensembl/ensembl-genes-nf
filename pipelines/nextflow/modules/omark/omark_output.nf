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


process OMARK_OUTPUT {
    tag "omark_output:${insdc_acc}"
    label 'local'
    publishDir "${params.outDir}/${publish_dir_name}/statistics", mode: 'copy'

    input:
        tuple val(insdc_acc), val(dbname), val(formated_sci_name), val(publish_dir_name), path(summary_file), val(omark_dir)

    output:
        tuple val(insdc_acc), val(dbname), path("*.txt")

    script:
        accession_formatted = insdc_acc.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")
        summary_name = [formated_sci_name.toLowerCase(), accession_formatted, "omark", "proteins_detailed_summary.txt"].join("_")

        """
        cat ${summary_file} > ${summary_name}
        """
}
