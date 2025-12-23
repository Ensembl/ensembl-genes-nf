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



process POPULATE_DB {
    label 'default'
    tag "$meta.core"

    input:
    //tuple val(gca), val(core), path(sql_file)
    tuple val(meta), path(sql_file)
    output:
    path("versions.yml"), emit: versions_file, optional :true

    when:
    params.apply_ensembl_stats || params.apply_beta_metakeys
    
    script:
    
    //${params.host} -w ${core} < ${statistics_file}
    // /hps/software/users/ensembl/ensw/mysql-cmds/ensembl/ensadmin/mysql-ens-genebuild-prod-6 ftricomi_gca035666275v1_core_110 </hps/nobackup/flicek/ensembl/genebuild/ftricomi/aves/chukar_partridge_annotation/alectoris_chukar/GCA_035666275.1//stats_ftricomi_gca035666275v1_core_110.sql
    """
    ${params.mysql_ensadmin}/${params.host} ${meta.core} < ${sql_file}
    """
}
