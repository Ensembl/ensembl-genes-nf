#!/usr/bin/env nextflow
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

include { getMetaValue } from './utils.nf'

process COPY_OUTPUT_TO_ENSEMBL_FTP {
    tag "copy on ftp"
    label 'default'

    input:
    tuple val(gca), val(dbname),path(summary_file)

    script:
    scientific_name = getMetaValue(dbname, "species.scientific_name")[0].meta_value.toString().replaceAll("\\s", "_")
    species=scientific_name.toLowerCase()
    gca_string = gca.toLowerCase().replaceAll(/\./, "v").replaceAll(/_/, "")
    publish_dir =scientific_name +'/'+gca+'/'+getMetaValue(dbname, "genebuild.annotation_source")[0].meta_value.toString() 
    statistics_files = "${params.outDir}/$publish_dir/statistics/*summary.txt"
    ftp_stats = "${params.production_ftp_dir}/$publish_dir/statistics" 
    ftp_path = "${params.production_ftp_dir}/$scientific_name"
    """
    sudo -u genebuild mkdir -p $ftp_stats; \
    sudo -u genebuild cp -f ${params.readme} $ftp_stats; 
    sudo -u genebuild cp -f $statistics_files  $ftp_stats; \
    sudo -u genebuild chmod 775 $ftp_stats/* -R;
    sudo -u genebuild chgrp ensemblftp $ftp_stats/* -R;
    """
    //sudo -u genebuild rsync -ahvW $summary_file $ftp_stats && rsync -avhc $summary_file $ftp_stats; \
}
