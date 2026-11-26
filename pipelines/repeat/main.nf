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
    ENSEMBL REPEAT ANNOTATION PIPELINE
========================================================================================

    This Nextflow DSL2 workflow performs repeat annotation of genome assemblies
    using RepeatModeler, RepeatMasker, RED, DustMasker, and Tandem Repeat Finder
    (TRF).

    The pipeline supports reuse of existing repeat annotations and repeat
    libraries, de novo RepeatModeler library generation, and optional publication
    of generated results to the Ensembl FTP infrastructure.

    Pipeline Stages:
    ----------------

    1. FETCH_GENOME
       Retrieve or locate genome assemblies specified in the input CSV.

    2. CHECK_REPEATMASKER_OUTPUT
       Check whether RepeatMasker annotations already exist on Ensembl FTP.
       Assemblies with existing annotations skip RepeatMasker processing.

    3. REPEAT LIBRARY SELECTION
       For assemblies requiring RepeatMasker annotation:
       - Use a custom library when supplied in the input CSV.
       - Otherwise, check for an existing RepeatModeler library.

    4. CHECK_REPEATMODELER_LIBRARY
       Check whether a precomputed RepeatModeler library is available.

    5a. DOWNLOAD_REPEATMODELER_LIBRARY
        Download an existing RepeatModeler library when available
        and RepeatMasker processing is enabled.

    5b. GENERATE_REPEATMODELER_LIBRARY
        Generate a de novo RepeatModeler library when no existing
        library is available and generate_lib is enabled.

    6. UPLOAD_REPEAT_LIBRARY_INTO_FTP
       Upload newly generated RepeatModeler libraries to Ensembl FTP.

    7. RUN_REPEATMASKER
       Identify and mask repetitive sequences using the selected library.
       Runs only for assemblies missing existing RepeatMasker annotations.

    8. RUN_RED
       Independently identify repetitive regions using RED.

    9. RUN_DUST
       Independently identify low-complexity regions using DustMasker.

    10. RUN_TRF
        Independently identify tandem repeats using TRF.

    11. UPLOAD_REPEATS_INTO_FTP
        Optionally upload newly generated RepeatMasker outputs
        to the configured Ensembl FTP destination.

    12. COLLECT_SOFTWARE_VERSIONS
        Collect and consolidate software version information
        from the executed modules.

    Input:
    ------

    CSV metadata file containing:
        - gca
        - species_name
        - genome_file
        - repeatmasker_library

    The genome_file and repeatmasker_library fields provide
    assembly and custom-library information where applicable.

    Output:
    -------

    Depending on the enabled analyses:
        - Retrieved genome assemblies
        - Downloaded or newly generated RepeatModeler libraries
        - RepeatMasker soft-masked genome outputs
        - RED repeat detection results
        - DustMasker low-complexity results
        - TRF tandem repeat results
        - Uploaded RepeatModeler libraries
        - Uploaded RepeatMasker annotations
        - Consolidated software version information

    Execution Parameters:
    ---------------------

        --csvFile          Input CSV metadata file
        --outdir           Pipeline output directory
        --generate_lib     Generate missing RepeatModeler libraries
        --run_repeatmasker Enable RepeatMasker analysis
        --run_red          Enable RED analysis
        --run_dust         Enable DustMasker analysis
        --run_trf          Enable TRF analysis
        --upload_repeats   Upload generated RepeatMasker outputs

    Workflow Notes:
    ---------------

    - Existing RepeatMasker annotations are checked before library selection.
    - Existing RepeatModeler libraries are reused when available.
    - Missing RepeatModeler libraries require generate_lib=true.
    - Custom RepeatMasker libraries bypass RepeatModeler library selection.
    - RED, DUST, and TRF operate independently on fetched genomes.
    - Only RepeatMasker outputs are currently connected to the repeat
        annotation upload channel.
    - Software versions are collected from executed processes.

========================================================================================
*/

nextflow.enable.dsl=2
include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'


// Import process modules
include { FETCH_GENOME }                     from './modules/fetch_genome.nf'
include { CHECK_FTP  as CHECK_REPEATMODELER_LIBRARY } from './modules/check_ftp.nf'
include { CHECK_FTP  as CHECK_REPEATMASKER_OUTPUT } from './modules/check_ftp.nf'
include { GENERATE_REPEATMODELER_LIBRARY }   from './modules/generate_repeatmodeler_library.nf'
include { DOWNLOAD_REPEATMODELER_LIBRARY }     from './modules/download_repeatmodeler_library.nf'
include { RUN_REPEATMASKER }                 from './modules/run_repeatmasker.nf'
include { RUN_RED }                          from './modules/run_red.nf'
include { RUN_DUST }                         from './modules/run_dust.nf'
include { RUN_TRF }                          from './modules/run_trf.nf'
include { UPLOAD_INTO_FTP as UPLOAD_REPEAT_LIBRARY_INTO_FTP }   from './modules/upload_into_ftp.nf'
include { UPLOAD_INTO_FTP as UPLOAD_REPEATS_INTO_FTP }          from './modules/upload_into_ftp.nf'
include { COLLECT_SOFTWARE_VERSIONS }        from './modules/collect_software_versions.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/




workflow REPEAT_ANNOTATION {

    take:
    csv_file

    main:

    // ========================================================
    // Initialize
    // ========================================================

    ch_versions_file = channel.empty()
    ch_repeat_output = channel.empty()

    // ========================================================
    // Stage 1: Read CSV
    // ========================================================

    data = channel.fromPath(
        csv_file,
        type: 'file',
        checkIfExists: true
    )
    .splitCsv(sep: ',', header: true)
    .map { row ->

        [
            gca: row.get('gca'),
            species_name: row.get('species_name'),
            genome_file: row.get('genome_file'),
            repeatmasker_library: row.get('repeatmasker_library')
        ]
    }

    // ========================================================
    // Stage 2: Fetch genomes
    // ========================================================

    genomeData = FETCH_GENOME(data).genome_file_output
        .map { meta, fa_file ->
            meta + [genome_file: fa_file]
        }

    ch_versions_file = ch_versions_file.mix(
        FETCH_GENOME.out.versions_file
    )
    // ========================================================
    // Stage 3: Check existing RepeatMasker outputs
    // ========================================================

    if (params.run_repeatmasker) {

        repeatmaskerCheckInput = genomeData.map { meta ->

            // IMPORTANT:
            // Replace this filename with the actual FTP output.
            def repeatmasker_url =
                "${params.repeatmasker_ftp_base}/${meta.species_name}/${meta.gca}_genome.masked.fa"

            tuple(meta, repeatmasker_url, "repeatmasker")
        }

        checkedRepeatMasker = CHECK_REPEATMASKER_OUTPUT(
            repeatmaskerCheckInput
        ).status_out

        ch_versions_file = ch_versions_file.mix(
            CHECK_REPEATMASKER_OUTPUT.out.versions_file
        )

        checkedRepeatMasker
            .map { meta, _url, _file_type, status_file ->
            def status = status_file.text.trim()
        def exists = status.equalsIgnoreCase("true")

        log.info "RepeatMasker check: ${meta.taxon_id} status=${status} exists=${exists}"

        tuple(meta, exists)
                //tuple(meta, status_file.text.trim() == "true")
            }
            .branch { _meta, exists ->
                available: exists
                missing: !exists
            }
            .set { repeatmasker_status }
        repeatmasker_status.available.view {
        meta, _exists ->
            "SKIPPING RepeatMasker: ${meta.taxon_id} (already available)"
        }

        repeatmasker_status.missing.view {
            meta, _exists ->
                "RUNNING RepeatMasker: ${meta.taxon_id} (missing)"
        }
        // Only assemblies missing RepeatMasker annotations
        // should proceed to library selection.
        repeatmaskerData = repeatmasker_status.missing
            .map { meta, _exists -> meta }

    } else {

        // No RepeatMasker checks needed.
        repeatmaskerData = channel.empty()
    }
    // ========================================================
    // Stage 4: Decide which libraries are needed
    // ========================================================

    // If RepeatMasker is enabled and the CSV provides a
    // library, use that library directly.
    //
    // Otherwise, use the RepeatModeler library workflow.

    if (params.run_repeatmasker) {

        repeatmaskerData
            .branch { meta ->

                custom_library:
                    meta.repeatmasker_library?.trim()

                needs_repeatmodeler:
                    !meta.repeatmasker_library?.trim()
            }
            .set { input_library_status }

        // ----------------------------------------------------
        // User-provided RepeatMasker libraries
        // ----------------------------------------------------

        customLibraries = input_library_status.custom_library
            .map { meta ->

                def library = file(
                    meta.repeatmasker_library,
                    checkIfExists: true
                )

                tuple(meta, library)
            }

        repeatModelerData = input_library_status.needs_repeatmodeler

    } else {

        customLibraries = channel.empty()

        // Generate/check RepeatModeler libraries independently
        // when RepeatMasker is disabled but generation requested.
        repeatModelerData = params.generate_lib ?
            genomeData : channel.empty()
    }

    // ========================================================
    // Stage 5: Check existing RepeatModeler libraries
    // ========================================================

    repeatModelInput = repeatModelerData
        .map { meta ->

            def repeat_url =
                "${params.repeatmodeler_ftp_base}/${meta.species_name}/${meta.gca}.repeatmodeler.fa"

            tuple(meta, repeat_url, "repeatmodeler")
        }

    checkedLibraries = CHECK_REPEATMODELER_LIBRARY(
        repeatModelInput
    ).status_out

    ch_versions_file = ch_versions_file.mix(
        CHECK_REPEATMODELER_LIBRARY.out.versions_file
    )

    // ========================================================
    // Stage 6: Branch based on remote library availability
    // ========================================================

    checkedLibraries
        .map { meta, repeat_url, _file_type, status_file ->

            def available = status_file.text.trim() == "true"

            tuple(meta, repeat_url, available)
        }
        .branch { _meta, _repeat_url, available ->

            available: available

            missing: !available
        }
        .set { library_status }

    // ========================================================
    // Stage 7a: Download existing RepeatModeler libraries
    // ========================================================
    // Download existing RepeatModeler libraries only if RepeatMasker is enabled
    rmlibrary = channel.empty()
    if (params.run_repeatmasker) {

        downloadInput = library_status.available
            .map { meta, repeat_url, _available ->

                tuple(repeat_url, meta)
            }

        rmlibrary = DOWNLOAD_REPEATMODELER_LIBRARY(
            downloadInput
        ).repeatmodeler_library_out

        ch_versions_file = ch_versions_file.mix(
            DOWNLOAD_REPEATMODELER_LIBRARY.out.versions_file
        )
    }
    // ========================================================
    // Stage 7b: Generate missing RepeatModeler libraries
    // ========================================================

    generateLibraryInput = library_status.missing
        .map { meta, _repeat_url, _available ->

            if (!params.generate_lib) {
                error (
                    "No RepeatModeler library available for " +
                    "${meta.gca}, and params.generate_lib is false"
                )
            }

            meta
        }

    ftpInput = GENERATE_REPEATMODELER_LIBRARY(
        generateLibraryInput
    ).repeatmodeler_library_out
    

    ch_versions_file = ch_versions_file.mix(
        GENERATE_REPEATMODELER_LIBRARY.out.versions_file
    )

    // ========================================================
    // Stage 7: Upload newly generated RepeatModeler libraries
    // ========================================================

    repeatModelerInput = UPLOAD_REPEAT_LIBRARY_INTO_FTP(
        ftpInput, params.repeatmodeler_ftp_local
    ).library_out

    ch_versions_file = ch_versions_file.mix(
        UPLOAD_REPEAT_LIBRARY_INTO_FTP.out.versions_file
    )

    // ========================================================
    // Stage 9: Merge RepeatModeler libraries
    // ========================================================

    allRepeatModelerLibraries = repeatModelerInput
        .mix(rmlibrary)

    // ========================================================
    // Stage 9: Select libraries for RepeatMasker
    // ========================================================

    if (params.run_repeatmasker) {

        // Priority:
        // 1. CSV repeatmasker_library
        // 2. Downloaded RepeatModeler library
        // 3. Newly generated RepeatModeler library

        allLibraries = customLibraries
            .mix(allRepeatModelerLibraries)

        RUN_REPEATMASKER(allLibraries)

        ch_versions_file = ch_versions_file.mix(
            RUN_REPEATMASKER.out.versions_file
        )
        repeatmasker_output = RUN_REPEATMASKER.out.repeatmasker_softmasked_out
        ch_repeat_output = ch_repeat_output.mix(
            repeatmasker_output
        )
    }

    // ========================================================
    // Stage 10: RED
    // ========================================================

    if (params.run_red) {

        RUN_RED(genomeData)

        ch_versions_file = ch_versions_file.mix(
            RUN_RED.out.versions_file
        )

        //ch_repeat_output = ch_repeat_output.mix(
        //    RUN_RED.out.repeat_output
        //)
    }

    // ========================================================
    // Stage 11: DUST
    // ========================================================

    if (params.run_dust) {

        RUN_DUST(genomeData)

        ch_versions_file = ch_versions_file.mix(
            RUN_DUST.out.versions_file
        )

        //ch_repeat_output = ch_repeat_output.mix(
        //   RUN_DUST.out.repeat_output
        //)
    }

    // ========================================================
    // Stage 12: TRF
    // ========================================================

    if (params.run_trf) {

        RUN_TRF(genomeData)

        ch_versions_file = ch_versions_file.mix(
            RUN_TRF.out.versions_file
        )

        //ch_repeat_output = ch_repeat_output.mix(
        //    RUN_TRF.out.repeat_output
        //)
    }

    // ========================================================
    // Stage 13: Upload repeat annotations
    // ========================================================

    if (params.upload_repeats) {

        UPLOAD_REPEATS_INTO_FTP(ch_repeat_output, params.repeatmasker_ftp_local)

        ch_versions_file = ch_versions_file.mix(
            UPLOAD_REPEATS_INTO_FTP.out.versions_file
        )
    }

    // ========================================================
    // Stage 14: Collect software versions
    // ========================================================

    COLLECT_SOFTWARE_VERSIONS(
        ch_versions_file.collect()
    )
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




