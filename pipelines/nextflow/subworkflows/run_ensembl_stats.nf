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

include { RUN_STATISTICS } from '../modules/ensembl_statistics/run_statistics.nf'
include { CREATE_STATS_JSON } from '../modules/ensembl_statistics/create_stats_json.nf'
include { ADD_STATS_ON_CORE  } from '../modules/ensembl_statistics/add_stats_on_core.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_ENSEMBL_STATS{
    take:                 
    db_meta

    main:
//    def db_meta1=db_meta
  //  db_meta1.flatten().view { d -> "GCA1: ${d.gca}, Core name: ${d.core}"}
        def statisticsFile = RUN_STATISTICS (db_meta.flatten())
        //def(statistics_output,json_file) = CREATE_STATS_JSON(statisticsFile)
        if(params.apply_stats){
        ADD_STATS_ON_CORE(statisticsFile)
        }

}


