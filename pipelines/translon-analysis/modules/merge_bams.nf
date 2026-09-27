include { MERGE_RIBO_BAMS } from './inputs/merge_bams.nf'

workflow MERGE_RIBO_INPUTS {
    take:
    transcriptome
    genome
    offsets
    merge_inputs
    merge_group

    main:
    if (merge_inputs) {
        tx_grouped = transcriptome
            .map { meta, bam, bai -> tuple([id: meta.merge_group ?: merge_group, merge_group: meta.merge_group ?: merge_group, bam_type: 'transcriptome'], bam, bai) }
            .groupTuple(by: 0)
            .map { meta, bams, bais -> tuple(meta, bams, bais) }
        gn_grouped = genome
            .map { meta, bam, bai -> tuple([id: meta.merge_group ?: merge_group, merge_group: meta.merge_group ?: merge_group, bam_type: 'genome'], bam, bai) }
            .groupTuple(by: 0)
            .map { meta, bams, bais -> tuple(meta, bams, bais) }
        MERGE_RIBO_BAMS(tx_grouped.mix(gn_grouped))
        tx_out = MERGE_RIBO_BAMS.out.merged.filter { meta, _bam, _bai -> meta.bam_type == 'transcriptome' }
        gn_out = MERGE_RIBO_BAMS.out.merged.filter { meta, _bam, _bai -> meta.bam_type == 'genome' }
        offsets_out = channel.empty()
    } else {
        tx_out = transcriptome
        gn_out = genome
        offsets_out = offsets
    }

    emit:
    transcriptome = tx_out
    genome = gn_out
    offsets = offsets_out
    manifests = merge_inputs ? MERGE_RIBO_BAMS.out.manifest : channel.empty()
}
