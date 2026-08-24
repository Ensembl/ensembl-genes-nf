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
// nextflow-lint-disable-file
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_BUSCO } from './subworkflows/run_busco.nf'
include { RUN_OMARK } from './subworkflows/run_omark.nf'
include { RUN_ENSEMBL_STATS } from './subworkflows/run_ensembl_stats.nf'
include { COLLECT_SOFTWARE_VERSIONS } from './modules/collect_software_versions.nf'
include { PEPSTATS } from './subworkflows/pepstats.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper function to clean cache directory
def cleanCacheDirectory() {
    if (params.cleanCache) {
        try {
            def cacheDir = file(params.cacheDir)
            if (cacheDir.exists() && cacheDir.isDirectory()) {
                cacheDir.listFiles().each { f ->
                    if (f.isDirectory()) {
                        f.deleteDir()
                    } else {
                        f.delete()
                    }
                }
                log.info("Cleaning process completed successfully.")
            }
        } catch (Exception e) {
            log.error("Exception occurred while executing cleaning command: ${e.message}")
        }
    }
}
workflow {
    log.info("Pipeline started at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}")

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info(paramsSummaryLog(workflow))

    if (params.run_pepstats) {
        PEPSTATS(params.csvFile)
    }

    if (params.run_busco_core || params.run_busco_ncbi) {
        def share_proteins = params.run_pepstats && params.run_busco_core &&
            (params.busco_mode == 'protein' || params.busco_mode == 'both')
        RUN_BUSCO(
            params.csvFile,
            share_proteins ? PEPSTATS.out.proteins : channel.empty(),
            share_proteins
        )
    }

    if (params.run_omark) {
        RUN_OMARK(params.csvFile)
    }

    if (params.run_ensembl_stats || params.run_ensembl_beta_metakeys) {
        RUN_ENSEMBL_STATS(params.csvFile)
    }

    // Collect version outputs
    def busco_versions = params.run_busco_core || params.run_busco_ncbi ? RUN_BUSCO.out.versions : channel.empty()
    def omark_versions = params.run_omark ? RUN_OMARK.out.versions : channel.empty()
    def stats_versions = params.run_ensembl_stats || params.run_ensembl_beta_metakeys ? RUN_ENSEMBL_STATS.out.versions : channel.empty()
    def pepstats_versions = params.run_pepstats ? PEPSTATS.out.versions : channel.empty()

    // Mix all versions
    def ch_all_versions = channel.empty()
        .mix(busco_versions)
        .mix(omark_versions)
        .mix(stats_versions)
        .mix(pepstats_versions)

    // Merge into single file and publish
    COLLECT_SOFTWARE_VERSIONS(ch_all_versions.collect())
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
