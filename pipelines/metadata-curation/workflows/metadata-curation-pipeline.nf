#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXTRACT_SRA_METADATA     } from '../modules/data-source-extraction'
include { MAP_ONTOLOGY_TERMS       } from '../modules/metadata-standardization'
include { EXTRACT_STUDY_CONTEXT    } from '../modules/metadata-standardization'
include { PROCESS_SAMPLE_METADATA  } from '../modules/metadata-standardization'
include { DETECT_BIOGROUPS         } from '../modules/output-generation'
include { VALIDATE_QUALITY         } from '../modules/output-generation'
include { FORMAT_SAMPLE_TABLE      } from '../modules/output-generation'

workflow METADATA_CURATION_SIMPLE {

    main:

    ch_versions = Channel.empty()

    //
    // Extract metadata from SRA
    //
    EXTRACT_SRA_METADATA (
        params.query,
        params.max_results ?: 10,
        params.email ?: ""
    )
    ch_versions = ch_versions.mix(EXTRACT_SRA_METADATA.out.versions)

    // Create metadata channel with meta map
    ch_sra_meta = EXTRACT_SRA_METADATA.out.metadata.map { file -> 
        [ [id: 'sra', source: 'sra'], file ] 
    }

    //
    // Map terms to ontologies
    //
    ch_cache_dir = Channel.of(file("${params.outdir}/ontology_cache"))
    
    MAP_ONTOLOGY_TERMS (
        ch_sra_meta,
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

    //
    // Detect biogroups
    //
    DETECT_BIOGROUPS (
        PROCESS_SAMPLE_METADATA.out.processed_metadata
    )
    ch_versions = ch_versions.mix(DETECT_BIOGROUPS.out.versions)

    //
    // Validate quality
    //
    VALIDATE_QUALITY (
        PROCESS_SAMPLE_METADATA.out.processed_metadata,
        DETECT_BIOGROUPS.out.biogroups.map { meta, file -> file }
    )
    ch_versions = ch_versions.mix(VALIDATE_QUALITY.out.versions)

    //
    // Format sample table
    //
    FORMAT_SAMPLE_TABLE (
        PROCESS_SAMPLE_METADATA.out.processed_metadata,
        DETECT_BIOGROUPS.out.biogroups.map { meta, file -> file },
        params.output_format ?: "csv",
        false
    )
    ch_versions = ch_versions.mix(FORMAT_SAMPLE_TABLE.out.versions)

    emit:
    sample_table      = FORMAT_SAMPLE_TABLE.out.table
    biogroups         = DETECT_BIOGROUPS.out.biogroups
    validation_report = VALIDATE_QUALITY.out.report
    versions          = ch_versions
}

workflow {
    METADATA_CURATION_SIMPLE()
}