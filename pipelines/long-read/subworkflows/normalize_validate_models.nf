include { GTF_TO_BED12; GTF_TO_BED12 as ISOQUANT_TO_BED12; CANONICALISE_COMBINED_MODELS; VALIDATE_LONG_READ_MODELS; AUDIT_NATIVE_MODELS } from '../modules/normalization.nf'

workflow NORMALIZE_VALIDATE_MODELS {
    take:
    merged_models
    reference_dictionary

    main:
    gtf_inputs = merged_models.map { meta, backend, scope, gtf, filelist, filelist_checksum -> tuple(meta, backend, scope, gtf) }
    GTF_TO_BED12(gtf_inputs.filter { _meta, backend, _scope, _gtf -> backend != 'isoquant' })
    ISOQUANT_TO_BED12(gtf_inputs.filter { _meta, backend, _scope, _gtf -> backend == 'isoquant' })
    combined_bed12 = GTF_TO_BED12.out.bed12.mix(ISOQUANT_TO_BED12.out.bed12)
    CANONICALISE_COMBINED_MODELS(combined_bed12.map { meta, backend, scope, bed12, manifest -> tuple(meta, backend, scope, bed12, manifest) })
    VALIDATE_LONG_READ_MODELS(CANONICALISE_COMBINED_MODELS.out.canonical.map { meta, backend, scope, bed12, manifest -> tuple(meta, backend, scope, bed12, manifest) }, reference_dictionary)
    AUDIT_NATIVE_MODELS(VALIDATE_LONG_READ_MODELS.out.validated)

    emit:
    validated = VALIDATE_LONG_READ_MODELS.out.validated
    reports = GTF_TO_BED12.out.report.mix(ISOQUANT_TO_BED12.out.report)
    summary = AUDIT_NATIVE_MODELS.out.summary
    versions = channel.empty()
}
