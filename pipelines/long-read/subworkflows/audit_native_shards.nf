include { AUDIT_NATIVE_SHARD_MODELS; AUDIT_BACKEND_SHARD_PLAN; VALIDATE_BACKEND_SHARD_PLAN; EMIT_EMPTY_NATIVE_MODEL } from '../modules/normalization.nf'

workflow AUDIT_NATIVE_SHARDS {
    take:
    native_models
    empty_native_inputs
    shard_plans
    backends

    main:
    EMIT_EMPTY_NATIVE_MODEL(empty_native_inputs)
    all_native_models = native_models.mix(EMIT_EMPTY_NATIVE_MODEL.out.native_models.map { meta, backend_meta, product, status -> tuple(meta, backend_meta, product, status) })
    AUDIT_NATIVE_SHARD_MODELS(all_native_models)
    audited = AUDIT_NATIVE_SHARD_MODELS.out.audited
    native_audits = audited.map { _meta, _backend_meta, _product, _status, audit -> audit }.collect()
    AUDIT_BACKEND_SHARD_PLAN(shard_plans.map { meta, _bam, _bai, plan -> tuple(meta, plan) }, backends, native_audits)
    VALIDATE_BACKEND_SHARD_PLAN(AUDIT_BACKEND_SHARD_PLAN.out.audit)

    emit:
    audited = audited
    versions = channel.empty()
}
