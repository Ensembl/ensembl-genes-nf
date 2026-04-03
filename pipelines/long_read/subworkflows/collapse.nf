/*
 * COLLAPSE SUBWORKFLOW
 * Collapse overlapping long-read alignments into non-redundant transcript models.
 * All samples are processed independently; read support is pooled per locus
 * in the classify step.
 */

include { COLLAPSE_TRANSCRIPTS } from '../modules/collapse_transcripts.nf'

workflow COLLAPSE {
    take:
    bam_bai_ch      // [ meta, bam, bai ] — one entry per sample
    min_overlap     // val
    max_intron_size // val

    main:
    ch_versions = Channel.empty()

    // Collapse per-sample. For species with many samples this could be
    // extended to merge BAMs first; for now one GFF3 per sample is emitted
    // and the classify step aggregates evidence.
    COLLAPSE_TRANSCRIPTS(
        bam_bai_ch,
        min_overlap,
        max_intron_size
    )
    ch_versions = ch_versions.mix(COLLAPSE_TRANSCRIPTS.out.versions)

    emit:
    gff3     = COLLAPSE_TRANSCRIPTS.out.gff3   // [ meta, gff3 ]
    stats    = COLLAPSE_TRANSCRIPTS.out.stats  // [ meta, tsv  ]
    versions = ch_versions
}
