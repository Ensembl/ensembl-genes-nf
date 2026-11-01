include { STRINGTIE2_COLLAPSE; STRINGTIE3_COLLAPSE } from '../modules/models.nf'

workflow COLLAPSE_STRINGTIE {
    take:
    shards
    reference
    backends

    main:
    stringtie2_native_models = channel.empty()
    stringtie3_native_models = channel.empty()
    stringtie2_versions = channel.empty()
    stringtie3_versions = channel.empty()
    if (backends.contains('stringtie2')) {
        STRINGTIE2_COLLAPSE(shards, reference)
        stringtie2_native_models = STRINGTIE2_COLLAPSE.out.native_models
        stringtie2_versions = STRINGTIE2_COLLAPSE.out.versions
    }
    if (backends.contains('stringtie3')) {
        STRINGTIE3_COLLAPSE(shards, reference)
        stringtie3_native_models = STRINGTIE3_COLLAPSE.out.native_models
        stringtie3_versions = STRINGTIE3_COLLAPSE.out.versions
    }

    emit:
    native_models = stringtie2_native_models.mix(stringtie3_native_models).map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = stringtie2_versions.mix(stringtie3_versions)
}
