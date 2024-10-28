// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

// http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


nextflow.enable.dsl=2

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// VALIDATE INPUTS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

// Validate input parameters
validateParameters()
log.info paramsSummaryLog(workflow)

if (!params.enscode) {
    exit 1, "Undefined --enscode parameter. Please provide the enscode path"
}
if (!params.outDir) {
    exit 1, "Undefined --outDir parameter. Please provide the output directory's path"
}
if (!params.server_set){
    params.mysql_ensadmin = "ensadmin"
}
else{
    params.mysql_ensadmin = params.server_set
}

// if (params.csvFile) {
//     csvFile = file(params.csvFile, checkIfExists: true)
// } else {
//     exit 1, 'CSV file not specified!'
// }

busco_mode = []
if (params.busco_mode instanceof java.lang.String) {
    busco_mode = [params.busco_mode]
}
else {
    busco_mode = params.busco_mode
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_COREDB_METADATA } from '../modules/genome_metadata/prepare_metadata.nf'
include { RUN_BUSCO } from '../subworkflows/run_busco.nf'
include { RUN_OMARK } from '../subworkflows/run_omark.nf'
include { RUN_ENSEMBL_STATS } from '../subworkflows/run_ensembl_stats.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STATISTICS{
    if(params.run_busco_ncbi && !params.run_busco_core){
        // Read data from the CSV file, split it, and map each row to extract GCA and taxon values
        data = Channel.fromPath(params.csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> [gca:row.get('gca'), taxon_id:row.get('taxon_id'), core:'core']}
        def busco_mode = 'genome'
        def copyToFtp = false
        data1=data
        data1.each{ d-> d.view()}
        
        // Now run BUSCO on genme mode without database as input
        RUN_BUSCO(data, busco_mode, copyToFtp)
    }

    if (params.run_busco_core || params.run_omark || params.run_ensembl_stats) {
        
        
        rows = file(params.csvFile).readLines().drop(1) // Skip header
        rows.each { row ->
            def core = row
            
            genome_metadata = PREPARE_COREDB_METADATA(core, params.metatable_keys).core_metadata

            if (params.run_busco_core) {
                RUN_BUSCO(genome_metadata, busco_mode, params.copyToFtp)
            }

            if (params.run_omark) {
            RUN_OMARK(genome_metadata)
            }

            if (params.run_ensembl_stats || params.run_ensembl_beta_metakeys) {
            RUN_ENSEMBL_STATS(genome_metadata)
            }
        }
    }    
    
    if (params.cleanCache) {
        // Clean cache directories
        exec "rm -rf ${params.cacheDir}/*"
    }
    
}