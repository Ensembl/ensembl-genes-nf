nextflow.enable.dsl = 2

include { CANONICALISE_COMBINED_MODELS } from '../modules/canonicalise_models.nf'
include { VALIDATE_LONG_READ_MODELS } from '../modules/validate_models.nf'

workflow VALIDATE_COMBINED_MODELS {
    take:
    merged_bed

    main:
    // merged_bed: tuple val(backend), val(cohort_id), path(merged_models.bed)
    canonical_input = merged_bed.map { backend, cohort_id, bed -> tuple(backend, [id: cohort_id, backend: backend], bed) }
    CANONICALISE_COMBINED_MODELS(canonical_input)
    VALIDATE_LONG_READ_MODELS(CANONICALISE_COMBINED_MODELS.out.bed.map { backend, meta, bed -> tuple(backend, meta, bed) })

    emit:
    bed = CANONICALISE_COMBINED_MODELS.out.bed
    checksum = CANONICALISE_COMBINED_MODELS.out.checksum
    validation = VALIDATE_LONG_READ_MODELS.out.report
    versions = CANONICALISE_COMBINED_MODELS.out.versions.mix(VALIDATE_LONG_READ_MODELS.out.versions)
}
