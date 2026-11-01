include { SAMTOOLS_VIEW_SHARD; TAMA_PREPARE_SAM; TAMA_COLLAPSE; TAMA_BED12_TO_GTF; CLASSIFY_NATIVE_GTF; WRITE_MODEL_INPUT_CHECKSUMS } from '../modules/models.nf'

workflow COLLAPSE_TAMA {
    take:
    shards
    reference
    cap_mode
    end_mode
    run_mode

    main:
    SAMTOOLS_VIEW_SHARD(shards)
    tama_inputs = SAMTOOLS_VIEW_SHARD.out.sam.join(shards, by: [0, 1]).map { meta, shard_meta, alignment_sam, shard_bam, shard_bai -> tuple(meta, shard_meta, shard_bam, shard_bai, alignment_sam) }
    TAMA_PREPARE_SAM(tama_inputs)
    TAMA_COLLAPSE(TAMA_PREPARE_SAM.out.prepared, reference, cap_mode, end_mode, run_mode)
    TAMA_BED12_TO_GTF(TAMA_COLLAPSE.out.raw.map { meta, shard_meta, _bam, _bai, _sam, tama_bed -> tuple(meta, shard_meta, tama_bed) })
    CLASSIFY_NATIVE_GTF(TAMA_BED12_TO_GTF.out.gtf.map { meta, shard_meta, native_gtf -> tuple(meta, shard_meta, 'tama', native_gtf) }, 'TAMA produced no models')
    WRITE_MODEL_INPUT_CHECKSUMS(TAMA_COLLAPSE.out.raw.map { meta, shard_meta, bam, bai, sam, _bed -> tuple(meta, shard_meta, 'tama', [bam, bai, sam]) })

    emit:
    native_models = CLASSIFY_NATIVE_GTF.out.native_models.map { meta, shard_meta, backend, product, status -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: shard_meta.shard_id, native_format: 'gtf', model_count: 0, status: status.text.contains('EMPTY') ? 'EMPTY' : 'SUCCESS'], product, status) }
    versions = SAMTOOLS_VIEW_SHARD.out.versions.mix(TAMA_PREPARE_SAM.out.versions).mix(TAMA_COLLAPSE.out.versions).mix(TAMA_BED12_TO_GTF.out.versions).mix(CLASSIFY_NATIVE_GTF.out.versions).mix(WRITE_MODEL_INPUT_CHECKSUMS.out.versions)
    checksums = WRITE_MODEL_INPUT_CHECKSUMS.out.checksums
}
