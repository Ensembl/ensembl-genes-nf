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

nextflow.enable.dsl=2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_STATISTICS } from '../modules/run_statistics.nf'
include { RUN_ENSEMBL_META as RUN_BETA_METAKEYS } from '../modules/run_ensembl_meta.nf'
include { POPULATE_DB as ADD_STATS_ON_CORE  } from '../modules/populate_db.nf'
include { POPULATE_DB as ADD_BETA_UPDATES_ON_CORE  } from '../modules/populate_db.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_ENSEMBL_STATS{
    take:                 
    csvFile

    main:
    // Validate required parameters
    if (params.fetch == null) {
        error "params.fetch must be defined as true or false"
    }
    // Read data from the CSV file, split it, and map each row to extract GCA and taxon values
        def data = Channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> 
                    [gca:row.get('gca'), 
                    core:row.get('core'),
                    species_id:row.get('species_id')
                    ]
                    }
        if(params.run_ensembl_stats){
        def statisticsFile = RUN_STATISTICS(data).output.statistics_output
        ADD_STATS_ON_CORE(statisticsFile)

        }    
        def db_meta=data
        if(params.run_ensembl_beta_metakeys){
        def betaMetakeys = RUN_BETA_METAKEYS (db_meta)
        ADD_BETA_UPDATES_ON_CORE(betaMetakeys)
            }
        }


