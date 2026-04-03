/*
 * ALIGN_AND_FILTER SUBWORKFLOW
 * Run exonerate (cdna or protein model) on sequence batches,
 * filter hits by coverage/pid, emit per-batch GFF3.
 */

include { EXONERATE_CDNA    } from '../modules/exonerate_cdna.nf'
include { EXONERATE_PROTEIN } from '../modules/exonerate_protein.nf'
include { FILTER_EXONERATE  } from '../modules/filter_exonerate.nf'

workflow ALIGN_AND_FILTER {
    take:
    seq_batches   // [ meta, fasta ] — sequence batches (cdna or protein)
    genome_fasta  // path: softmasked genome
    query_type    // val: 'cdna' or 'protein'
    query_fasta   // path or null: original query FASTA for accurate coverage

    main:
    ch_versions = Channel.empty()

    if (query_type == 'cdna') {
        EXONERATE_CDNA(seq_batches, genome_fasta)
        ch_raw_gff = EXONERATE_CDNA.out.gff
        ch_versions = ch_versions.mix(EXONERATE_CDNA.out.versions)
    } else {
        EXONERATE_PROTEIN(seq_batches, genome_fasta)
        ch_raw_gff = EXONERATE_PROTEIN.out.gff
        ch_versions = ch_versions.mix(EXONERATE_PROTEIN.out.versions)
    }

    FILTER_EXONERATE(
        ch_raw_gff,
        query_type,
        query_fasta ?: []
    )
    ch_versions = ch_versions.mix(FILTER_EXONERATE.out.versions)

    emit:
    gff3     = FILTER_EXONERATE.out.gff3  // [ meta, filtered.gff3 ]
    versions = ch_versions
}
