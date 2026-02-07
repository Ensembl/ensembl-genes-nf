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
/*POPULATE_DB process to populate core database with statistics and metakeys
Inputs:
- meta: metadata map containing dbname, species_id, taxon_id, gca
- sql_file: path to the SQL file to be executed
Outputs:
- versions.yml: versions file containing software versions used
*/
process POPULATE_DB {
    label 'default'
    tag "${meta.dbname}"

    input:
    tuple val(meta), path(sql_file)

    output:
    path ("versions.yml"), emit: versions_file, optional: true

    when:
    params.apply_ensembl_stats || params.apply_ensembl_beta_metakeys

    script:
    """
    ${params.mysql_ensadmin}/${params.host} ${meta.dbname} < ${sql_file}

    # Create versions file
    MYSQL_VERSION=\$(mysql --version | grep -oP 'Distrib \\K[0-9.]+')
    
    echo '"POPULATE_DB":' > versions.yml
    echo "  mysql: \$MYSQL_VERSION" >> versions.yml
    """
}
