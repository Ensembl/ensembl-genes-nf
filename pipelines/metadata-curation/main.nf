#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Multi-Source Metadata Curation Pipeline
========================================================================================
    Github: https://github.com/ensembl-genes-nf
----------------------------------------------------------------------------------------
*/

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
if (!params.query) {
    error "Error: Please provide a search query using --query"
}

if (!params.outdir) {
    error "Error: Please provide an output directory using --outdir"
}


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { MULTI_SOURCE_INGESTION } from './subworkflows/multi-source-ingestion'
include { ONTOLOGY_STANDARDIZATION } from './subworkflows/ontology-standardization'  
include { CROSS_MODAL_ANALYSIS } from './subworkflows/cross-modal-analysis'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow METADATA_CURATION {

    take:
    query // channel: search query string

    main:

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Multi-source data ingestion
    //
    MULTI_SOURCE_INGESTION (
        query,
        params.max_results,
        params.include_pmc,
        params.email ?: ""
    )
    ch_versions = ch_versions.mix(MULTI_SOURCE_INGESTION.out.versions)

    //
    // SUBWORKFLOW: Ontology standardization
    //
    ONTOLOGY_STANDARDIZATION (
        MULTI_SOURCE_INGESTION.out.sra_metadata,
        MULTI_SOURCE_INGESTION.out.geo_metadata,
        MULTI_SOURCE_INGESTION.out.pmc_metadata,
        params.ontology_cache_dir ?: "${launchDir}/${params.outdir}/ontology_cache"
    )
    ch_versions = ch_versions.mix(ONTOLOGY_STANDARDIZATION.out.versions)
    ch_fixed_metadata = ONTOLOGY_STANDARDIZATION.out.processed_metadata
        .flatMap { items ->
            // Convert flattened list back to tuples
            def result = []
            for (int i = 0; i < items.size(); i += 2) {
                if (i + 1 < items.size()) {
                    result.add([items[i], items[i + 1]])
                }
            }
            return result
        }
    //
    // SUBWORKFLOW: Cross-modal analysis and biogroup detection
    //
    CROSS_MODAL_ANALYSIS (
        ch_fixed_metadata,
        params.output_format ?: "csv",
        params.output_level ?: "sample",
        params.include_raw_metadata ?: false
    )
    ch_versions = ch_versions.mix(CROSS_MODAL_ANALYSIS.out.versions)


    emit:
    sample_table      = CROSS_MODAL_ANALYSIS.out.sample_table
    biogroup_table    = CROSS_MODAL_ANALYSIS.out.biogroup_table
    validation_report = CROSS_MODAL_ANALYSIS.out.validation_report
    versions          = ch_versions
}

/*
========================================================================================
    THE END
========================================================================================
*/

workflow {
    query_ch = Channel.of(params.query)
    METADATA_CURATION(query_ch)
}