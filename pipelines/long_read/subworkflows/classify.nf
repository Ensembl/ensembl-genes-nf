/*
 * CLASSIFY SUBWORKFLOW
 * Run blastp against the protein database, then assign biotypes based on
 * protein support level.
 *
 * classify_transcripts.py handles ORF extraction from the GFF3 internally,
 * keeping the module interface simple. If ORF extraction is later separated
 * (e.g. via TransDecoder), split into: TRANSDECODER → BLAST_BLASTP → CLASSIFY.
 */

include { BLAST_BLASTP         } from '../modules/blast_blastp.nf'
include { CLASSIFY_TRANSCRIPTS } from '../modules/classify_transcripts.nf'

workflow CLASSIFY {
    take:
    gff3_ch    // [ meta, gff3 ]
    protein_db // [ meta, db_dir ]

    main:
    ch_versions = Channel.empty()

    BLAST_BLASTP(
        gff3_ch,
        protein_db,
        'tsv'
    )
    ch_versions = ch_versions.mix(BLAST_BLASTP.out.versions)

    CLASSIFY_TRANSCRIPTS(
        gff3_ch,
        BLAST_BLASTP.out.tsv
    )
    ch_versions = ch_versions.mix(CLASSIFY_TRANSCRIPTS.out.versions)

    emit:
    gff3     = CLASSIFY_TRANSCRIPTS.out.gff3  // [ meta, classified.gff3 ]
    versions = ch_versions
}
