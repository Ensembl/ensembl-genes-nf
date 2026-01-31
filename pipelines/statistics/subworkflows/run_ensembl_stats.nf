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

include { DB_METADATA } from '../modules/db_metadata.nf'
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
    // Read data from the CSV file, split it, and map each row to extract GCA and taxon values
        def data = Channel.fromPath(csvFile, type: 'file', checkIfExists: true)
                .splitCsv(sep:',', header:true)
                .map { row -> 
                    [gca:'UNKNOWN',
                    taxon_id:'UNKNOWN',
                    dbname:row.get('dbname'),
                    species_id:row.get('species_id')?row.get('species_id'):1
                    ]
                    }
        ch_versions_file = Channel.empty()            
        def metadata = DB_METADATA(data).metadata
        .map { meta, metadata_file ->
            def lines = metadata_file.text.readLines()
            def new_taxon = lines[0].split('=')[1]
            def new_gca = lines[1].split('=')[1]
            def production_name = lines[2].split('=')[1]
            def updated_meta = meta + [taxon_id: new_taxon, gca: new_gca, production_name: production_name]
            updated_meta
        }
        if(params.run_ensembl_stats){
        def statisticsFile = RUN_STATISTICS(metadata).statistics_output
        ch_versions_file = ch_versions_file.mix(RUN_STATISTICS.out.versions_file)
        ADD_STATS_ON_CORE(statisticsFile)
        ch_versions_file = ch_versions_file.mix(ADD_STATS_ON_CORE.out.versions_file)

        }    
        if(params.run_ensembl_beta_metakeys){
        def betaMetakeys = RUN_BETA_METAKEYS(metadata).ensembl_meta_output
        ch_versions_file = ch_versions_file.mix(RUN_BETA_METAKEYS.out.versions_file)
        ch_versions_file.view { "After RUN_BETA_METAKEYS mix: $it" }
        ADD_BETA_UPDATES_ON_CORE(betaMetakeys)
        ch_versions_file = ch_versions_file.mix(ADD_BETA_UPDATES_ON_CORE.out.versions_file)
        ch_versions_file.view { "After ADD_BETA_UPDATES_ON_CORE mix: $it" }
            }
            // Collect all versions to ensure they're ready before emitting
ch_versions_file = ch_versions_file.collect()
    emit:
    versions = ch_versions_file

        }


