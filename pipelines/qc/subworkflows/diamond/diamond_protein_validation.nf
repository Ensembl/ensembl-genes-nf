include { PREP_DIAMOND_DB } from '../../modules/diamond/prep_diamond_db.nf'
include { DIAMOND_BLASTP } from '../../modules/diamond/blastp.nf'
include { DIAMOND_PARSE } from '../../modules/diamond/parse_diamond.nf'

workflow DIAMOND_PROTEIN_VALIDATION {
    take:
        protein_ch
        annotation_ch
        ref_protein_faa
        ref_diamond_db

    main:
        if (ref_diamond_db) {
            diamond_db = channel.value(file(ref_diamond_db))
            database_versions = channel.empty()
        } else {
            reference_db = PREP_DIAMOND_DB(file(ref_protein_faa))
            diamond_db = reference_db.db.collect()
            database_versions = reference_db.versions
        }

        hits = DIAMOND_BLASTP(
            protein_ch,
            diamond_db
        )

        annotation_by_id = annotation_ch.map { meta, genome_fasta, gff3 ->
            tuple(meta.id, meta, genome_fasta, gff3)
        }
        hits_by_id = hits.hits.map { meta, diamond_hits ->
            tuple(meta.id, meta, diamond_hits)
        }
        parse_input = annotation_by_id.join(hits_by_id).map {
            sample_id, meta, genome_fasta, gff3, _hit_meta, diamond_hits ->
                tuple(meta, genome_fasta, gff3, diamond_hits)
        }
        parsed = DIAMOND_PARSE(parse_input)

    emit:
        diamond_tsv = hits.hits
        diamond_stats = parsed.stats
        versions = database_versions.mix(hits.versions).mix(parsed.versions)
}
