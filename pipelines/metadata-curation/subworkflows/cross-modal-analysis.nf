include { DETECT_BIOGROUPS      } from '../modules/detect_biogroups'
include { VALIDATE_QUALITY      } from '../modules/validate_quality'
include { FORMAT_SAMPLE_TABLE   } from '../modules/format_sample_table'
include { FORMAT_BIOGROUP_TABLE } from '../modules/format_biogroup_table'

workflow CROSS_MODAL_ANALYSIS {

    take:
    processed_metadata  // channel: [ [meta, processed_metadata.json], ... ]
    output_format       // val: output format (csv/tsv/json)
    output_level        // val: output level (sample/biogroup/study)
    include_raw         // val: include raw metadata and evidence details

    main:

    ch_versions = Channel.empty()
    //
    // Detect biogroups across all data
    //
    DETECT_BIOGROUPS (
        processed_metadata
    )
    ch_versions = ch_versions.mix(DETECT_BIOGROUPS.out.versions)

    //
    // Validate metadata quality
    //
    VALIDATE_QUALITY (
        processed_metadata,
        DETECT_BIOGROUPS.out.biogroups.map { meta, file -> file }
    )
    ch_versions = ch_versions.mix(VALIDATE_QUALITY.out.versions)

    //
    // Format sample-level output table
    //
    FORMAT_SAMPLE_TABLE (
        processed_metadata,
        DETECT_BIOGROUPS.out.biogroups.map { meta, file -> file },
        output_format,
        include_raw
    )
    ch_sample_table = FORMAT_SAMPLE_TABLE.out.table
    ch_versions = ch_versions.mix(FORMAT_SAMPLE_TABLE.out.versions)

    //
    // Format biogroup-level output table
    //
    FORMAT_BIOGROUP_TABLE (
        processed_metadata,
        DETECT_BIOGROUPS.out.biogroups.map { meta, file -> file },
        output_format,
        include_raw
    )
    ch_biogroup_table = FORMAT_BIOGROUP_TABLE.out.table
    ch_versions = ch_versions.mix(FORMAT_BIOGROUP_TABLE.out.versions)

    emit:
    sample_table      = ch_sample_table                    // path: output_table.*
    biogroup_table    = ch_biogroup_table                  // path: output_table.*  
    biogroups         = DETECT_BIOGROUPS.out.biogroups    // channel: [ [meta, biogroups.json] ]
    biogroup_report   = DETECT_BIOGROUPS.out.report       // path: biogroup_report.json
    validation_report = VALIDATE_QUALITY.out.report       // path: validation_report.json
    versions          = ch_versions                        // channel: [ versions.yml ]
}