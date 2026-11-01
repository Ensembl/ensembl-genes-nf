include { SAMTOOLS_VIEW_SHARD; BAM_TO_ALIGNMENT_GTF; TMERGE_COLLAPSE } from '../modules/models.nf'

workflow COLLAPSE_TMERGE {
    take:
    shards
    tmerge_args

    main:
    SAMTOOLS_VIEW_SHARD(shards)
    BAM_TO_ALIGNMENT_GTF(SAMTOOLS_VIEW_SHARD.out.sam)
    TMERGE_COLLAPSE(BAM_TO_ALIGNMENT_GTF.out.alignment_gtf, tmerge_args)

    emit:
    native_models = TMERGE_COLLAPSE.out.native_models.map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = SAMTOOLS_VIEW_SHARD.out.versions.mix(BAM_TO_ALIGNMENT_GTF.out.versions).mix(TMERGE_COLLAPSE.out.versions)
}
