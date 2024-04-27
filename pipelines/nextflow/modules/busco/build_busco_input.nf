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
@Grab('org.codehaus.groovy:groovy-all:2.2.2')
import groovy.sql.Sql
include { getMetaValue } from '../utils.nf'

process BUILD_BUSCO_INPUT {

    label 'default'
    tag "prepare input data for $core"

    input:
    tuple val(gca), val(taxon_id), val(core)

    output:
    tuple val(out_gca), val(out_taxon_id), val(core)

    script:
    println(core)
    def gca_value = gca
    def taxon_id_value = taxon_id

    if (gca=='gca'){
       // query = "SELECT meta_value FROM meta WHERE meta_key= 'assembly.accession'"
         meta_gca = getMetaValue(core, "assembly.accession")
         //gca1.view()
         //gca =channel.sql.fromQuery(query, db: 'core_db')
         gca_value=meta_gca[0].meta_value.toString()
    //     out_gca = "'${gca_value}'" 
//         println(out_gca)
    }
    
    if (taxon_id=='taxon_id'){
        meta_taxon_id = getMetaValue(core, "species.taxonomy_id") 
        taxon_id_value = meta_taxon_id[0].meta_value.toString()
      //  out_taxon_id = "'${taxon_id_value}'"

    }
    out_gca = gca_value
    out_taxon_id = taxon_id_value
}
