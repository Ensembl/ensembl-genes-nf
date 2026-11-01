include { ISOQUANT_ANNOTATION_FREE; ISOQUANT_REFERENCE_GUIDED; NORMALISE_ISOQUANT_OUTPUTS; AUDIT_ISOQUANT_OUTPUTS } from '../modules/optional_backends.nf'

workflow RUN_ISOQUANT {
    take:
    alignment
    reference
    mode
    genedb
    data_type
    analysis
    high_memory
    isoquant_args
    large_output

    main:
    native_models = channel.empty()
    versions = channel.empty()
    if (mode == 'annotation_free') {
        ISOQUANT_ANNOTATION_FREE(alignment, reference, data_type, analysis, high_memory, isoquant_args, large_output)
        NORMALISE_ISOQUANT_OUTPUTS(ISOQUANT_ANNOTATION_FREE.out.native_models)
        AUDIT_ISOQUANT_OUTPUTS(NORMALISE_ISOQUANT_OUTPUTS.out.native_models)
        native_models = AUDIT_ISOQUANT_OUTPUTS.out.native_models
        versions = ISOQUANT_ANNOTATION_FREE.out.versions
    } else {
        ISOQUANT_REFERENCE_GUIDED(alignment, reference, genedb, data_type, analysis, high_memory, isoquant_args, large_output)
        NORMALISE_ISOQUANT_OUTPUTS(ISOQUANT_REFERENCE_GUIDED.out.native_models)
        AUDIT_ISOQUANT_OUTPUTS(NORMALISE_ISOQUANT_OUTPUTS.out.native_models)
        native_models = AUDIT_ISOQUANT_OUTPUTS.out.native_models
        versions = ISOQUANT_REFERENCE_GUIDED.out.versions
    }
    emit:
    native_models = native_models.map { meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: 'merged', native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = versions.mix(NORMALISE_ISOQUANT_OUTPUTS.out.versions).mix(AUDIT_ISOQUANT_OUTPUTS.out.versions)
}
