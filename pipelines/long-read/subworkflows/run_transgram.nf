include { TRANSGRAM_ANNOTATION_FREE; TRANSGRAM_ANNOTATION_GUIDED; FINALIZE_TRANSGRAM_GTF } from '../modules/optional_backends.nf'

workflow RUN_TRANSGRAM {
    take:
    alignment
    annotation
    data_type
    transgram_args
    transgram_executable

    main:
    typed_alignment = alignment.map { meta, bam, bai -> tuple(meta, bam, bai, meta.transgram_data_type ?: data_type) }
    native_models = channel.empty()
    versions = channel.empty()
    if (annotation) {
        TRANSGRAM_ANNOTATION_GUIDED(typed_alignment, annotation, transgram_args, transgram_executable)
        FINALIZE_TRANSGRAM_GTF(TRANSGRAM_ANNOTATION_GUIDED.out.raw)
        native_models = FINALIZE_TRANSGRAM_GTF.out.native_models
        versions = TRANSGRAM_ANNOTATION_GUIDED.out.versions
    } else {
        TRANSGRAM_ANNOTATION_FREE(typed_alignment, transgram_args, transgram_executable)
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
