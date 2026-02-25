/*
 * QUALITY CONTROL SUBWORKFLOW
 * Validates collapsed reads using getRPF
 */

include { CHECK_CLEANLINESS } from '../modules/check_cleanliness.nf'

workflow QUALITY_CONTROL {
    take:
    collapsed_reads   // tuple: [ meta, collapsed_fasta ]

    main:
    // Check cleanliness of collapsed reads
    CHECK_CLEANLINESS(
        collapsed_reads,
        'RPF'  // count_pattern parameter
    )

    // Pass through the collapsed reads for downstream processing
    clean_samples = collapsed_reads

    emit:
    samples = clean_samples                    // tuple: [ meta, collapsed_fasta ]
    reports = CHECK_CLEANLINESS.out.report     // tuple: [ meta, report ]
    rpf_checks = CHECK_CLEANLINESS.out.rpf_checks  // tuple: [ meta, checks ]
}
