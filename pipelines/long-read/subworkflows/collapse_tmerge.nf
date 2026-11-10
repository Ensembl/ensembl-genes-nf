include { SAMTOOLS_VIEW_SHARD; BAM_TO_ALIGNMENT_GTF; TMERGE_COLLAPSE; CLASSIFY_NATIVE_GTF; WRITE_MODEL_INPUT_CHECKSUMS } from '../modules/models.nf'

workflow COLLAPSE_TMERGE {
    take:
    shards
    tmerge_args

    main:
    SAMTOOLS_VIEW_SHARD(shards)
    BAM_TO_ALIGNMENT_GTF(SAMTOOLS_VIEW_SHARD.out.sam)
    TMERGE_COLLAPSE(BAM_TO_ALIGNMENT_GTF.out.alignment_gtf, tmerge_args)
    CLASSIFY_NATIVE_GTF(TMERGE_COLLAPSE.out.raw.map { meta, shard_meta, _alignment_gtf, native_gtf -> tuple(meta, shard_meta, 'tmerge', native_gtf) }, 'TMerge produced no models')
    WRITE_MODEL_INPUT_CHECKSUMS(TMERGE_COLLAPSE.out.raw.map { meta, shard_meta, alignment_gtf, _native_gtf -> tuple(meta, shard_meta, 'tmerge', [alignment_gtf]) })

    emit:
    native_models = CLASSIFY_NATIVE_GTF.out.native_models.map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = SAMTOOLS_VIEW_SHARD.out.versions.mix(TMERGE_COLLAPSE.out.versions)
}
