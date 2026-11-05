include { TRANSGRAM_ANNOTATION_FREE; TRANSGRAM_ANNOTATION_GUIDED; FINALIZE_TRANSGRAM_GTF } from '../modules/optional_backends.nf'

workflow RUN_TRANSGRAM {
    take:
    alignment
    annotation
    data_type
    transgram_args
    transgram_executable

    main:
    native_models = channel.empty()
    versions = channel.empty()
    if (annotation) {
        TRANSGRAM_ANNOTATION_GUIDED(alignment, annotation, data_type, transgram_args, transgram_executable)
        FINALIZE_TRANSGRAM_GTF(TRANSGRAM_ANNOTATION_GUIDED.out.raw)
        native_models = FINALIZE_TRANSGRAM_GTF.out.native_models
        versions = TRANSGRAM_ANNOTATION_GUIDED.out.versions
    } else {
        TRANSGRAM_ANNOTATION_FREE(alignment, data_type, transgram_args, transgram_executable)
        FINALIZE_TRANSGRAM_GTF(TRANSGRAM_ANNOTATION_FREE.out.raw)
        native_models = FINALIZE_TRANSGRAM_GTF.out.native_models
        versions = TRANSGRAM_ANNOTATION_FREE.out.versions
    }

    emit:
    native_models = native_models.map { meta, backend, product, status ->
        tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: 'merged', native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status)
    }
    versions = versions
}
