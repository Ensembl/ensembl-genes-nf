include { STRINGTIE2_COLLAPSE; STRINGTIE3_COLLAPSE; CLASSIFY_NATIVE_GTF as CLASSIFY_STRINGTIE2; CLASSIFY_NATIVE_GTF as CLASSIFY_STRINGTIE3; WRITE_MODEL_INPUT_CHECKSUMS as CHECKSUM_STRINGTIE2; WRITE_MODEL_INPUT_CHECKSUMS as CHECKSUM_STRINGTIE3 } from '../modules/models.nf'

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
        CLASSIFY_STRINGTIE2(STRINGTIE2_COLLAPSE.out.raw.map { meta, shard_meta, _bam, _bai, native_gtf -> tuple(meta, shard_meta, 'stringtie2', native_gtf) }, 'StringTie2 produced no models')
        CHECKSUM_STRINGTIE2(STRINGTIE2_COLLAPSE.out.raw.map { meta, shard_meta, bam, bai, _gtf -> tuple(meta, shard_meta, 'stringtie2', [bam, bai]) })
        stringtie2_native_models = CLASSIFY_STRINGTIE2.out.native_models
        stringtie2_versions = STRINGTIE2_COLLAPSE.out.versions.mix(CLASSIFY_STRINGTIE2.out.versions).mix(CHECKSUM_STRINGTIE2.out.versions)
    }
    if (backends.contains('stringtie3')) {
        STRINGTIE3_COLLAPSE(shards, reference)
        CLASSIFY_STRINGTIE3(STRINGTIE3_COLLAPSE.out.raw.map { meta, shard_meta, _bam, _bai, native_gtf -> tuple(meta, shard_meta, 'stringtie3', native_gtf) }, 'StringTie3 produced no models')
        CHECKSUM_STRINGTIE3(STRINGTIE3_COLLAPSE.out.raw.map { meta, shard_meta, bam, bai, _gtf -> tuple(meta, shard_meta, 'stringtie3', [bam, bai]) })
        stringtie3_native_models = CLASSIFY_STRINGTIE3.out.native_models
        stringtie3_versions = STRINGTIE3_COLLAPSE.out.versions.mix(CLASSIFY_STRINGTIE3.out.versions).mix(CHECKSUM_STRINGTIE3.out.versions)
    }

    emit:
    native_models = stringtie2_native_models.mix(stringtie3_native_models).map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = stringtie2_versions.mix(stringtie3_versions)
}
