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
    tag "Load_stats:${dbname}"

    input:
        tuple val(insdc_acc), val(taxonomy_id), val(dbname), 
            val(production_name), val(organism_name), val(annotation_source), path(sql_file)

    script:
        if (params.server_set){
            host_admin = params.server_set
        }
        else{
            host_admin = params.mysql_ensadmin
        }
        """
        ${params.mysql_cmds}/${host_admin}/${params.host} ${dbname} < ${sql_file}
        """
}