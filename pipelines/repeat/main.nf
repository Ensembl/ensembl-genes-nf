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
    4. RUN_REPEATMASKER       - Identify and mask repeats using combined libraries
    
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

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

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
                    genome_file: row.get('genome_file')
                ]
            }
    
    
        // Stage 1: Fetch genome assemblies from NCBI
        //FETCH_GENOME(data)
        genomeData = FETCH_GENOME(data).genome_file_output
            .map { meta, fna_file ->
                return tuple(meta, fna_file)
            }
            .view { meta, fna_file -> "Processing: gca=${meta.gca},  genome=${fna_file}" }
        ch_versions_file = ch_versions_file.mix(FETCH_GENOME.out.versions_file)
        if (params.generate_lib) {

        // Stage 2: Check for existing RepeatModeler libraries
        checkedLibraries = FETCH_REPEAT_MODEL(genomeData).rep_library_file_output
            .view { meta, genome_file, library_file -> "Checked for RepeatModeler library for ${meta.gca}, file: ${library_file}" }
        ch_versions_file = ch_versions_file.mix(FETCH_REPEAT_MODEL.out.versions_file)

        // Separate genomes based on library availability
        checkedLibraries
            .branch { _meta, _genome_file, library_file ->
                // Check if library file contains error message
                available: !library_file.text.contains("No repeatmodeler file available")
                missing: library_file.text.contains("No repeatmodeler file available")
            }
            .set { library_status }
        
        // Stage 3a: Generate de novo RepeatModeler libraries for genomes without existing libraries
        repeatModelerInput = library_status.missing
            .map { meta, genome_file, _library_file -> tuple(meta, genome_file) }
        GENERATE_REPEATMODELER_LIBRARY(repeatModelerInput)
        ch_versions_file = ch_versions_file.mix(GENERATE_REPEATMODELER_LIBRARY.out.versions_file)
        
        
        // Stage 3b: Download pre-computed libraries for genomes that have them
        // Transform to [url, gca] format expected by CHECK_AND_DOWNLOAD_RMLIBRARY
        downloadInput = library_status.available
            .map { meta, genome_file, _library_file ->
                def url = "${params.repeats_ftp_base}/${meta.species_name}/${meta.gca}.repeatmodeler.fa"
                return tuple(url, meta, genome_file)
            }
        CHECK_AND_DOWNLOAD_RMLIBRARY(downloadInput)
        ch_versions_file = ch_versions_file.mix(CHECK_AND_DOWNLOAD_RMLIBRARY.out.versions_file)
                
        // Merge both library sources (generated + downloaded)
        // Outputs have format: tuple val(meta), path(genome_file), path(library_file)
        allLibraries = GENERATE_REPEATMODELER_LIBRARY.out.repeatmodeler_library_out
            .mix(CHECK_AND_DOWNLOAD_RMLIBRARY.out.repeatmodeler_library_out)
            .view { meta, genome_file, library -> "Library ready for ${meta.gca}: ${library}" }
        
        if(params.run_repeatmasker) {
        
        // Stage 4: Run RepeatMasker to identify and mask repeats
        RUN_REPEATMASKER(allLibraries)
        
        ch_versions_file = ch_versions_file.mix(RUN_REPEATMASKER.out.versions_file)
        }
        }
        if (params.run_red) {
            // Run RED for repeat annotation
            // Similar approach: join genome files with RED results
            RUN_RED(genomeData)
            ch_versions_file = ch_versions_file.mix(RUN_RED.out.versions_file)
        }
        if (params.run_dust) {
            // Run DUST for repeat annotation
            // Similar approach: join genome files with DUST results
            RUN_DUST(genomeData) 
            ch_versions_file = ch_versions_file.mix(RUN_DUST.out.versions_file)          
        }
        if (params.run_trf) {
            // Run TRF for repeat annotation
            // Similar approach: join genome files with TRF results
            RUN_TRF(genomeData)
            ch_versions_file = ch_versions_file.mix(RUN_TRF.out.versions_file)  
        }

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
    // Workflow execution handlers
    //workflow.onStart {
    //    log.info """
    //    ================================================================================
    //    REPEAT ANNOTATION PIPELINE
     //   ================================================================================
    //    Output directory : ${params.outdir}
     //   CSV input file   : ${params.csvFile ?: 'NOT PROVIDED'}
    //    NCBI base URL    : ${params.ncbiBaseUrl}
    //    Repeats FTP base : ${params.repeats_ftp_base}
    //    RepeatModeler    : ${params.repeatmodeler_path}
     //   RepeatMasker     : ${params.repeatmasker_path}
     //   ================================================================================
    //    """.stripIndent()
        
        // Validate required parameters
    //    if (!params.outdir) {
    //        error "ERROR: --outdir parameter is required. Please provide the output directory path."
    //    }
        
    //    if (!params.csvFile) {
    //        error "❌ ERROR: --csvFile parameter is required. Please provide the path to the CSV file."
    //    }
        
        // Check if CSV file exists
    //    if (!file(params.csvFile).exists()) {
    //        error "❌ ERROR: CSV file does not exist: ${params.csvFile}"
    //    }
        
    //    log.info "✅ Parameters validated successfully"
    //}
    
    workflow.onComplete {
        log.info """
        ================================================================================
        Pipeline Execution Summary
        ================================================================================
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        Exit status  : ${workflow.exitStatus}
        Work directory: ${workflow.workDir}
        ================================================================================
        """.stripIndent()
    }
    
    workflow.onError {
        log.error """
        ================================================================================
        Pipeline Execution Error
        ================================================================================
        Error message: ${workflow.errorMessage}
        Error report : ${workflow.errorReport}
        ================================================================================
        """.stripIndent()
    }
    

    // Execute main workflow
    REPEAT_ANNOTATION(params.csvFile)
}
