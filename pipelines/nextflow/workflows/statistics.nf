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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


if (!params.enscode) {
    exit 1, "Undefined --enscode parameter. Please provide the enscode path"
}
if (!params.outDir) {
    exit 1, "Undefined --outDir parameter. Please provide the output directory's path"
}

if (params.csvFile) {
    csvFile = file(params.csvFile, checkIfExists: true)
} else {
    exit 1, 'CSV file not specified!'
}

busco_mode = []
if (params.busco_mode instanceof java.lang.String) {
    busco_mode = [params.busco_mode]
}
else {
    busco_mode = params.busco_mode
}

acceptable_projects = ['ensembl', 'brc']
if (!acceptable_projects.contains(params.project)) {
    exit 1, 'Invalid project name'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info ''
    log.info 'Pipeline to run Busco score in protein and/or genome mode'
    log.info '-------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow -C ensembl-genes-nf/pipelines/nextflow/workflows/nextflow.config \
                run ensembl-genes-nf/pipelines/nextflow/workflows/main.nf \
                -entry STATISTICS --enscode --csvFile --outDir --host --port --user --bioperl --project \
                --run_busco_core --run_busco_ncbi --run_omark --run_ensembl_stats \
                --apply_stats --copyToFtp --busco_mode --apply_busco_metakeys'
    log.info ''
    log.info 'Options:'
    log.info '  --host STR                   Db host server'
    log.info '  --port INT                   Db port'
    log.info '  --user STR                   Db user'
    log.info '  --user_r STR                 Db user read_only'
    log.info '  --enscode STR                Enscode path'
    log.info '  --outDir STR                 Output directory.'
    log.info '  --csvFile STR                Path for the csv containing the db name' 
    log.info '  --bioperl STR                BioPerl path (optional)'
    log.info '  --project STR                Project, for the formatting of the output ("ensembl" or "brc")'
    log.info '  --run_busco_core bool        Run BUSCO (protein or genome mode see --busco_mode) given a mysql db, default false'
    log.info '  --run_busco_ncbi bool        Run BUSCO given a assembly_accession and taxonomy id in genome mode only, default false'
    log.info '  --run_omark bool             Run OMARK given a mysql db, default false'
    log.info '  --run_ensembl_stats bool     Run Ensembl statistics given a mysql db, default false'
    log.info '  --apply_stats bool           Upload Ensembl statistics in a mysql db, default false'
    log.info '  --copyToFtp bool             Copy output in Ensembl ftp, default false'
    log.info '  --busco_mode STR             Busco mode: genome or protein, default is to run both'
    log.info '  --busco_dataset STR          Busco dataset: optional'
    log.info '  --apply_busco_metakeys bool  Create JSON file with Busco metakeys and load it into the db, default false'

    exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_BUSCO } from '../subworkflows/run_busco.nf'
include { RUN_OMARK } from '../subworkflows/run_omark.nf'
include { RUN_ENSEMBL_STATS } from '../subworkflows/run_ensembl_stats.nf'

include { getMetaValue } from '../modules/utils.nf'
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
        RUN_BUSCO(data, busco_mode, copyToFtp)
    }
    if (params.run_busco_core || params.run_omark || params.run_ensembl_stats) {
    
        def fullProcessedData = []
        rows = file(params.csvFile).readLines().drop(1) // Skip header
        rows.each { row ->
            def core = row
            def gcaChannel = getMetaValue(core, "assembly.accession")
            def taxonIdChannel = getMetaValue(core, "species.taxonomy_id")
            def gca = gcaChannel[0].meta_value.toString()
            def taxon_id = taxonIdChannel[0].meta_value.toString()
            fullProcessedData.add([gca: gca, taxon_id: taxon_id, core: core])
        }

    full_data = Channel.from(fullProcessedData)
    
    def minimalProcessedData = []
        rows = file(params.csvFile).readLines().drop(1) // Skip header
        rows.each { row ->
            def core = row
            def gcaChannel = getMetaValue(core, "assembly.accession")
            def gca = gcaChannel[0].meta_value.toString()
            minimalProcessedData.add([gca: gca, core: core])
        }

    minimal_data = Channel.from(minimalProcessedData)

        if (params.run_busco_core) {
        RUN_BUSCO(full_data, busco_mode, params.copyToFtp)
        }

        if (params.run_omark) {
        RUN_OMARK(minimal_data)
        }

        if (params.run_ensembl_stats) {
        RUN_ENSEMBL_STATS(minimal_data)
        }
    }    
    
    if (params.cleanCache) {
        // Clean cache directories
        exec "rm -rf ${params.cacheDir}/*"
    }
    
}
