include { PREP_DIAMOND_DB } from '../../modules/diamond/prep_diamond_db.nf'
include { DIAMOND_BLASTP } from '../../modules/diamond/blastp.nf'

workflow DIAMOND_PROTEIN_VALIDATION {
    take:
        protein_ch
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

    emit:
        diamond_tsv = hits.hits
        versions = database_versions.mix(hits.versions)
}
