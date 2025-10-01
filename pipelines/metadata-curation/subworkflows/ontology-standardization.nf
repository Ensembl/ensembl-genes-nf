include { MAP_ONTOLOGY_TERMS    } from '../modules/map_ontology_terms'
include { EXTRACT_STUDY_CONTEXT } from '../modules/extract_study_context'
include { PROCESS_SAMPLE_METADATA } from '../modules/process_sample_metadata'

workflow ONTOLOGY_STANDARDIZATION {

    take:
    sra_metadata        // path: sra_metadata.json
    geo_metadata        // path: geo_metadata.json  
    pmc_metadata        // path: pmc_content.json
    ontology_cache_dir  // val: directory for ontology cache

    main:

    ch_versions = Channel.empty()

    // Create metadata channel with meta map
    ch_sra_meta = sra_metadata.map { file -> 
        [ [id: 'sra', source: 'sra'], file ] 
    }
    ch_geo_meta = geo_metadata.map { file -> 
        [ [id: 'geo', source: 'geo'], file ] 
    }
    
    // Combine all metadata sources
    ch_all_metadata = ch_sra_meta.mix(ch_geo_meta)

    //
    // Create ontology cache directory as string parameter
    //
    ch_cache_dir = Channel.value(ontology_cache_dir)
    
    //
    // Map terms to ontologies
    // 
    MAP_ONTOLOGY_TERMS (
        ch_all_metadata,
        ch_cache_dir
    )
    ch_versions = ch_versions.mix(MAP_ONTOLOGY_TERMS.out.versions)

    //
    // Extract study-level context
    //
    EXTRACT_STUDY_CONTEXT (
        MAP_ONTOLOGY_TERMS.out.mapped_metadata
    )
    ch_versions = ch_versions.mix(EXTRACT_STUDY_CONTEXT.out.versions)

    //
    // Process sample-level metadata with context
    //
    PROCESS_SAMPLE_METADATA (
        EXTRACT_STUDY_CONTEXT.out.context_metadata
    )
    ch_versions = ch_versions.mix(PROCESS_SAMPLE_METADATA.out.versions)

    emit:
    processed_metadata = PROCESS_SAMPLE_METADATA.out.processed_metadata.collect() // channel: [ [meta, processed_metadata.json], ... ]
    versions          = ch_versions                                               // channel: [ versions.yml ]
}