/*
 * SEARCH_NCRNA SUBWORKFLOW
 * Run cmsearch (Rfam) and optional BLASTN (miRBase) on genome chunks,
 * filter hits, and emit per-chunk GFF3.
 */

include { CMSEARCH    } from '../modules/cmsearch.nf'
include { BLASTN_MIRNA } from '../modules/blastn_mirna.nf'
include { FILTER_NCRNA } from '../modules/filter_ncrna.nf'

workflow SEARCH_NCRNA {
    take:
    genome_chunks   // [ meta, fasta ] — one per genomic chunk
    rfam_cm         // path: Rfam.cm covariance model file
    mirna_fasta     // path or null: miRBase all_mirnas.fa (set null to skip)
    genome_blast_db // [ meta, db_files ] or null: BLAST DB of genome

    main:
    ch_versions = Channel.empty()

    // Attach CM file to each chunk for cmsearch
    ch_cmsearch_in = genome_chunks.map { meta, fasta -> [meta, rfam_cm, fasta] }

    CMSEARCH(ch_cmsearch_in)
    ch_versions = ch_versions.mix(CMSEARCH.out.versions)

    // Optional miRNA BLAST
    ch_blast_tsv = Channel.empty()
    if (mirna_fasta && genome_blast_db) {
        ch_mirna = Channel.of([[id: 'mirna'], file(mirna_fasta)])
        BLASTN_MIRNA(ch_mirna, genome_blast_db)
        ch_versions = ch_versions.mix(BLASTN_MIRNA.out.versions)
        ch_blast_tsv = BLASTN_MIRNA.out.tsv.map { meta, tsv -> tsv }.first()
    }

    // Filter per chunk (pass blast TSV or [] if not available)
    FILTER_NCRNA(
        CMSEARCH.out.tblout,
        ch_blast_tsv.ifEmpty([])
    )
    ch_versions = ch_versions.mix(FILTER_NCRNA.out.versions)

    emit:
    gff3     = FILTER_NCRNA.out.gff3  // [ meta, ncrna.gff3 ]
    versions = ch_versions
}
