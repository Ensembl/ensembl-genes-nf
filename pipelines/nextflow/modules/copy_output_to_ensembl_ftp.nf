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

process COPY_OUTPUT_TO_ENSEMBL_FTP {
    tag "${formated_sci_name}:${insdc_acc}"
    label 'local'

    input:
        tuple val(insdc_acc), val(dbname), val(formated_sci_name), 
            val(publish_dir_name), path(summary_file)

    script:
        statistics_files = file("${params.outDir}/${publish_dir_name}/statistics/*summary.txt")
        read_me = file("${workflow.projectDir}/../data/README.txt")
        ftp_stats = "${params.production_ftp_dir}/${publish_dir_name}/statistics"
        // ftp_path = "${params.production_ftp_dir}/$formated_sci_name"
        """
        sudo -u genebuild mkdir -p ${ftp_stats}; \
        sudo -u genebuild cp -f ${read_me} ${ftp_stats};
        sudo -u genebuild cp -f ${statistics_files} ${ftp_stats}; \
        sudo -u genebuild chmod 775 ${ftp_stats}/* -R;
        sudo -u genebuild chgrp ensemblftp ${ftp_stats}/* -R;
        """
        //sudo -u genebuild rsync -ahvW $summary_file $ftp_stats && rsync -avhc $summary_file $ftp_stats; \
}