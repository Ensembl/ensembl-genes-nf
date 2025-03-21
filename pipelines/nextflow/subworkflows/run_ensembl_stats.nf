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

include { POPULATE_DB as ADD_BETA_UPDATES_ON_CORE  } from '../modules/ensembl_statistics/populate_db.nf'
include { POPULATE_DB as ADD_STATS_ON_CORE  } from '../modules/ensembl_statistics/populate_db.nf'
include { RUN_ENSEMBL_META as RUN_BETA_METAKEYS } from '../modules/ensembl_statistics/run_ensembl_meta.nf'
include { RUN_STATISTICS } from '../modules/ensembl_statistics/run_statistics.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_ENSEMBL_STATS{
    take:                 
        db_meta

    main:
        if( params.run_ensembl_stats ) {

            statisticsFile = RUN_STATISTICS(db_meta)

            if( params.apply_ensembl_stats ) {
                ADD_STATS_ON_CORE(statisticsFile)
                }
            }

        if(params.run_ensembl_beta_metakeys){

            betaMetakeys = RUN_BETA_METAKEYS (db_meta)

            if(params.apply_ensembl_beta_metakeys){
                ADD_BETA_UPDATES_ON_CORE(betaMetakeys)
                }
        }
}


