include { BAMBU_DISCOVERY } from '../modules/optional_backends.nf'

workflow RUN_BAMBU {
    take:
    alignment
    reference
    track_reads
    quantify

    main:
    BAMBU_DISCOVERY(alignment, reference, track_reads, quantify)
    emit:
    native_models = BAMBU_DISCOVERY.out.native_models.map { meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: 'merged', native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = BAMBU_DISCOVERY.out.versions
}
