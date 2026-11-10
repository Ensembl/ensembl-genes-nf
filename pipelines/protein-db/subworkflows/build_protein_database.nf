include { FETCH_UNIPROT_FASTA } from '../modules/fetch_uniprot_fasta.nf'
include { BUILD_DIAMOND_DB } from '../modules/build_diamond_db.nf'

workflow BUILD_PROTEIN_DATABASE {
    take:
    builds_ch
    uniprot_api
    page_size
    fetch_script

    main:
    // Public boundary: [meta, build], where build contains source, scope,
    // taxon_id, query_extra, and include_isoforms.
    FETCH_UNIPROT_FASTA(builds_ch, uniprot_api, page_size, fetch_script)
    BUILD_DIAMOND_DB(FETCH_UNIPROT_FASTA.out.fasta)

    emit:
    databases = BUILD_DIAMOND_DB.out.database
    versions = FETCH_UNIPROT_FASTA.out.versions.mix(BUILD_DIAMOND_DB.out.versions)
}
