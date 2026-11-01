include { COMPARE_NORMALIZED_MODELS } from '../modules/normalization.nf'

workflow COMPARE_MODELS {
    take:
    records
    cohort_id

    main:
    COMPARE_NORMALIZED_MODELS(records, cohort_id)

    emit:
    comparison = COMPARE_NORMALIZED_MODELS.out.comparison
    versions = COMPARE_NORMALIZED_MODELS.out.versions
}
