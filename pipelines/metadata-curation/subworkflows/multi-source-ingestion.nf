include { EXTRACT_SRA_METADATA } from '../modules/extract_sra_metadata'
include { EXTRACT_GEO_METADATA } from '../modules/extract_geo_metadata'
include { EXTRACT_PMC_CONTENT  } from '../modules/extract_pmc_content'

workflow MULTI_SOURCE_INGESTION {

    take:
    query           // val: search query string
    max_results     // val: maximum results per source
    include_pmc     // val: whether to include PMC content
    email           // val: email for API requests

    main:

    ch_versions = Channel.empty()

    //
    // Extract metadata from SRA
    //
    EXTRACT_SRA_METADATA (
        query,
        max_results,
        email
    )
    ch_versions = ch_versions.mix(EXTRACT_SRA_METADATA.out.versions)

    //
    // Extract metadata from GEO  
    //
    EXTRACT_GEO_METADATA (
        query,
        max_results,
        email
    )
    ch_versions = ch_versions.mix(EXTRACT_GEO_METADATA.out.versions)

    //
    // Extract content from PMC (optional)
    //
    ch_pmc_metadata = Channel.empty()
    if (include_pmc) {
        EXTRACT_PMC_CONTENT (
            query,
            max_results,
            email,
            true  // include_full_text
        )
        ch_versions = ch_versions.mix(EXTRACT_PMC_CONTENT.out.versions)
        ch_pmc_metadata = EXTRACT_PMC_CONTENT.out.content
    }

    emit:
    sra_metadata = EXTRACT_SRA_METADATA.out.metadata     // path: sra_metadata.json
    geo_metadata = EXTRACT_GEO_METADATA.out.metadata     // path: geo_metadata.json
    pmc_metadata = ch_pmc_metadata                        // path: pmc_content.json
    versions     = ch_versions                            // channel: [ versions.yml ]
}