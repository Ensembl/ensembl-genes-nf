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

/*
========================================================================================
    REPEAT ANNOTATION PIPELINE
========================================================================================
    This workflow performs comprehensive repeat annotation on genome assemblies using
    RepeatModeler for library generation and RepeatMasker/DustMasker/TRF for repeat identification.
    

    Pipeline Stages:
    1. FETCH_GENOME           - Download genome assemblies from NCBI
    2. FETCH_REPEAT_MODEL     - Check for existing RepeatModeler libraries
    3a. GENERATE_REPEATMODELER_LIBRARY - Generate de novo libraries for genomes without existing ones
    3b. CHECK_AND_DOWNLOAD_RMLIBRARY   - Download pre-computed libraries when available
    3c. UPLOAD_REPEAT_LIBRARY_INTO_FTP - Upload generated libraries to Ensembl FTP server
    4. RUN_REPEATMASKER       - Identify and mask repeats using combined libraries
    5. RUN_RED                 - Identify repetitive regions using RED
    6. RUN_DUST                - Identify repetitive regions using DUST
    7. RUN_TRF                 - Identify repetitive regions using TRF
    8. UPLOAD_REPEATS_INTO_FTP - Upload RepeatModeler libraries and RepeatMasker results to Ensembl FTP server
    9. COLLECT_SOFTWARE_VERSIONS - Collect software versions from all processes and merge them into a single file
    Input:
        CSV file with columns: species_name, GCA_accession
    
    Output:
        - Downloaded genome assemblies
        - RepeatModeler libraries (generated or downloaded)
        - RepeatMasker annotation results (masked sequences, GTF files, statistics)
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'


// Import process modules
include { FETCH_GENOME }                     from './modules/fetch_genome.nf'
include { FETCH_REPEAT_MODEL }               from './modules/fetch_repeat_model.nf'
include { GENERATE_REPEATMODELER_LIBRARY }   from './modules/generate_repeatmodeler_library.nf'
include { CHECK_AND_DOWNLOAD_RMLIBRARY }     from './modules/check_and_download_rmlibrary.nf'
include { RUN_REPEATMASKER }                 from './modules/run_repeatmasker.nf'
include { RUN_RED }                          from './modules/run_red.nf'
include { RUN_DUST }                         from './modules/run_dust.nf'
include { RUN_TRF }                          from './modules/run_trf.nf'
include { UPLOAD_REPEAT_LIBRARY_INTO_FTP }   from './modules/upload_repeat_library_into_ftp.nf'
include { UPLOAD_REPEATS_INTO_FTP }          from './modules/upload_repeats_into_ftp.nf'
include { COLLECT_SOFTWARE_VERSIONS }        from './modules/collect_software_versions.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
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

workflow REPEAT_ANNOTATION {
    take:
    csv_file
    main:
    
        // Initialize versions channel
        ch_versions_file = channel.empty()
        data = channel.fromPath(csv_file, type: 'file', checkIfExists: true)
            .splitCsv(sep: ',', header: true)
            .map { row ->
                [
                    gca: row.get('gca'),
                    species_name: row.get('species_name'),
                    genome_file: row.get('genome_file'),
                    repeatmasker_library: row.get('repeatmasker_library')
                ]
            }
    
    
        // Stage 1: Fetch genome assemblies from NCBI
        //FETCH_GENOME(data)
        genomeData = FETCH_GENOME(data).genome_file_output
            .map { meta, fna_file ->
                return meta + [ genome_file: fna_file ]
            }
            .view { meta -> "Processing: gca=${meta.gca},${meta.species_name},  genome=${meta.genome_file}" }
        ch_versions_file = ch_versions_file.mix(FETCH_GENOME.out.versions_file)
        if (params.generate_lib) {

        // Stage 2: Check for existing RepeatModeler libraries
        checkedLibraries = FETCH_REPEAT_MODEL(genomeData).rep_library_file_output
            .view { meta, library_file -> "Checked for RepeatModeler library for ${meta.gca}, genome_file ${meta.genome_file}, file: ${library_file}" }
        ch_versions_file = ch_versions_file.mix(FETCH_REPEAT_MODEL.out.versions_file)

        // Separate genomes based on library availability
        checkedLibraries
        .view { meta,  library_file ->
                "meta=${meta}  library=${library_file}"
                    }
            .branch { _meta,  library_file ->
                // Check if library file contains error message
                available: !library_file.text.contains("No repeatmodeler file available")
                missing: library_file.text.contains("No repeatmodeler file available")
            }
            .set { library_status }
        
        // Stage 3a: Generate de novo RepeatModeler libraries for genomes without existing libraries
        generateLibraryInput = library_status.missing
            .map { meta,  _library_file -> meta }
        ftpInput=GENERATE_REPEATMODELER_LIBRARY(generateLibraryInput).repeatmodeler_library_out
        ch_versions_file = ch_versions_file.mix(GENERATE_REPEATMODELER_LIBRARY.out.versions_file)
        
        repeatModelerInput = UPLOAD_REPEAT_LIBRARY_INTO_FTP(ftpInput).library_out
        ch_versions_file = ch_versions_file.mix(UPLOAD_REPEAT_LIBRARY_INTO_FTP.out.versions_file) 
        // Stage 3b: Download pre-computed libraries for genomes that have them
        // Transform to [url, gca] format expected by CHECK_AND_DOWNLOAD_RMLIBRARY
        downloadInput = library_status.available
            .map { meta,  _library_file ->
                def url = "${params.repeats_ftp_base}/${meta.species_name}/${meta.gca}.repeatmodeler.fa"
                return tuple(url, meta)
            }
        // nextflow-lint-disable    
        rmlibrary = CHECK_AND_DOWNLOAD_RMLIBRARY(downloadInput).repeatmodeler_library_out // nextflow-lint-disable
        ch_versions_file = ch_versions_file.mix(CHECK_AND_DOWNLOAD_RMLIBRARY.out.versions_file)
        if(params.run_repeatmasker) {        
        // Merge both library sources (generated + downloaded)
        allLibraries = repeatModelerInput
            .mix(rmlibrary)
            .view { meta,  library -> "Library ready for ${meta.gca}: ${library}" }
        ch_repeat_output = channel.empty()
        
            // Stage 4: Run RepeatMasker to identify and mask repeats
            RUN_REPEATMASKER(allLibraries)
            ch_versions_file = ch_versions_file.mix(RUN_REPEATMASKER.out.versions_file)
            ch_repeat_output = ch_repeat_output.mix(RUN_REPEATMASKER.out.repeatmasker_out)
            }
        }
        if (params.run_red) {
            // Run RED for repeat annotation
            RUN_RED(genomeData)
            ch_versions_file = ch_versions_file.mix(RUN_RED.out.versions_file)
            ch_repeat_output = ch_repeat_output.mix(RUN_RED.out.repeat_output)

        }
        if (params.run_dust) {
            // Run DUST for repeat annotation
            RUN_DUST(genomeData) 
            ch_versions_file = ch_versions_file.mix(RUN_DUST.out.versions_file)   
            ch_repeat_output = ch_repeat_output.mix(RUN_DUST.out.repeat_output)       
        }
        if (params.run_trf) {
            // Run TRF for repeat annotation
            RUN_TRF(genomeData)
            ch_versions_file = ch_versions_file.mix(RUN_TRF.out.versions_file)  
            ch_repeat_output = ch_repeat_output.mix(RUN_TRF.out.repeat_output)

        }
        
        if (params.upload_repeats) {
        UPLOAD_REPEATS_INTO_FTP(ch_repeat_output)
        ch_versions_file = ch_versions_file.mix(UPLOAD_REPEATS_INTO_FTP.out.versions_file)
        }

        // Merge into single file and publish
        COLLECT_SOFTWARE_VERSIONS(ch_versions_file.collect())

}

/*
========================================================================================
    ENTRY WORKFLOW
========================================================================================
*/

workflow {
    log.info("Pipeline started at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}")
    // Validate input parameters
    validateParameters()
    // Print summary of supplied parameters
    log.info(paramsSummaryLog(workflow))
    // Execute main workflow
    REPEAT_ANNOTATION(params.csvFile)
}

//onComplete {
 //   log.info("Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}")
 //   log.info("Execution status: ${workflow.success ? 'Successful' : 'Failed'}")
 //   cleanCacheDirectory()
//}


