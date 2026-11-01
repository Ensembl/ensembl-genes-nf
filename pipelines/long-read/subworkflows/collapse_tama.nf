include { TAMA_PREPARE_SAM; TAMA_COLLAPSE } from '../modules/models.nf'

workflow COLLAPSE_TAMA {
    take:
    shards
    reference
    cap_mode
    end_mode
    run_mode

    main:
    TAMA_PREPARE_SAM(shards)
    TAMA_COLLAPSE(TAMA_PREPARE_SAM.out.prepared, reference, cap_mode, end_mode, run_mode)

    emit:
    native_models = TAMA_COLLAPSE.out.native_models.map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = TAMA_PREPARE_SAM.out.versions.mix(TAMA_COLLAPSE.out.versions)
    checksums = TAMA_COLLAPSE.out.checksums
}
