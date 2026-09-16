nextflow.enable.dsl = 2

include { CANONICALISE_COMBINED_MODELS } from '../modules/canonicalise_models.nf'
include { VALIDATE_LONG_READ_MODELS } from '../modules/validate_models.nf'

workflow VALIDATE_COMBINED_MODELS {
    take:
    merged_bed
    cohort_id

    main:
    // merged_bed: path merged TAMA BED
    canonical_input = merged_bed.map { bed -> tuple([id: cohort_id], bed) }
    CANONICALISE_COMBINED_MODELS(canonical_input)
    VALIDATE_LONG_READ_MODELS(CANONICALISE_COMBINED_MODELS.out.bed.map { _meta, bed -> bed })

    emit:
    bed = CANONICALISE_COMBINED_MODELS.out.bed
    checksum = CANONICALISE_COMBINED_MODELS.out.checksum
    validation = VALIDATE_LONG_READ_MODELS.out.report
    versions = CANONICALISE_COMBINED_MODELS.out.versions
}
