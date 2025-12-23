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

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters;paramsSummaryLog } from 'plugin/nf-schema'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_BUSCO } from 'subworkflows/run_busco.nf'
include { RUN_OMARK } from 'subworkflows/run_omark.nf'
include { RUN_ENSEMBL_STATS } from 'subworkflows/run_ensembl_stats.nf'

//include { getMetaValue } from '../modules/utils.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def deleteRecursively(Path path) {
    if (java.nio.file.Files.isDirectory(path)) {
        java.nio.file.Files.newDirectoryStream(path).each { subPath ->
            deleteRecursively(subPath)
        }
    }
    java.nio.file.Files.delete(path)
}

workflow {
    log.info "Pipeline started at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
    // Validate input parameters
    validateParameters()
    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)
            if (params.run_busco_core) {
        RUN_BUSCO(params.csvFile)
        }
        if (params.run_omark) {
        RUN_OMARK(params.csvFile)
        }
        if (params.run_ensembl_stats || params.run_ensembl_beta_metakeys) {
        RUN_ENSEMBL_STATS(params.csvFile)
        }



workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
    log.info  "Execution status: ${workflow.success ? 'Successful' : 'Failed'}"
    if (params.cleanCache) {    
    try {
        def outDir = java.nio.file.Paths.get(params.cacheDir)
        java.nio.file.Files.newDirectoryStream(outDir, "*").each { path ->
        if (!path.toString().endsWith(".gz")) {
            deleteRecursively(path)
            }
        }
        log.info "Cleaning process completed successfully."
    } catch (Exception e) {
        log.error "Exception occurred while executing cleaning command: ${e.message}", e
    }
    }
}

workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

}