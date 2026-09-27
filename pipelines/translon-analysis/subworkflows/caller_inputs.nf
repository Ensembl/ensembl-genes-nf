include { PREPARE_TRANSCRIPT_MODELS } from '../modules/preparation/transcript_models.nf'
include { PREPARE_RIBORF_READS } from '../modules/preparation/riborf_offset_correct.nf'
include { PREPARE_RIBOTRICER_OFFSETS } from '../modules/preparation/ribotricer_offsets.nf'
include { PREPARE_RIBOTRICER_AUTO_OFFSETS } from '../modules/preparation/ribotricer_auto_offsets.nf'

workflow PREPARE_CALLER_INPUTS {
    take:
    transcriptome
    gtf
    offsets

    main:
    PREPARE_TRANSCRIPT_MODELS(transcriptome, gtf)
    models_with_offsets = PREPARE_TRANSCRIPT_MODELS.out.models
        .map { meta, bam, bai, genepred, bed12 -> tuple(meta.id, meta, bam, bai, genepred, bed12) }
        .join(offsets.map { meta, offset -> tuple(meta.id, offset) }, by: 0)
        .map { _id, meta, bam, bai, genepred, bed12, offset -> tuple(meta, bam, bai, genepred, bed12, offset) }
    PREPARE_RIBORF_READS(models_with_offsets)

    emit:
    transcript_models = PREPARE_TRANSCRIPT_MODELS.out.models
    riborf = PREPARE_RIBORF_READS.out.riborf
}
