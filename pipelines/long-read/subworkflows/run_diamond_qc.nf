include { EXTRACT_COMBINED_TRANSCRIPTS; PREDICT_LONGEST_ATG_ORFS; PUBLISH_ORF_PRODUCTS; PREP_DIAMOND_DB; WRITE_DIAMOND_DATABASE_IDENTITY as WRITE_PROTEIN_IDENTITY; WRITE_DIAMOND_DATABASE_IDENTITY as WRITE_DATABASE_IDENTITY; DIAMOND_BLASTP; REPORT_COMBINED_DIAMOND_MODELS; COMPARE_CANDIDATE_MODELS } from '../modules/finalization.nf'
include { PREPARE_DIAMOND_QUERY; EMIT_EMPTY_DIAMOND_HITS } from '../modules/diamond.nf'

workflow RUN_DIAMOND_QC {
    take:
    validated_models
    reference_proteins
    diamond_database
    diamond_enabled
    orf_enabled

    main:
    EXTRACT_COMBINED_TRANSCRIPTS(validated_models)
    PREDICT_LONGEST_ATG_ORFS(EXTRACT_COMBINED_TRANSCRIPTS.out.transcripts)
    orf_products = channel.empty()
    orf_versions = channel.empty()
    if (orf_enabled) {
        PUBLISH_ORF_PRODUCTS(PREDICT_LONGEST_ATG_ORFS.out.orfs)
        orf_products = PUBLISH_ORF_PRODUCTS.out.orfs
        orf_versions = PUBLISH_ORF_PRODUCTS.out.versions
    }
    database_versions = channel.empty()
    diamond_versions = channel.empty()
    diamond_reports = channel.empty()
    comparisons = channel.empty()
    if (diamond_enabled && diamond_database) {
        database = diamond_database
    } else if (diamond_enabled) {
        PREP_DIAMOND_DB(reference_proteins)
        database = PREP_DIAMOND_DB.out.database
        WRITE_PROTEIN_IDENTITY(reference_proteins, 'reference_proteins')
        database_versions = PREP_DIAMOND_DB.out.versions.mix(WRITE_PROTEIN_IDENTITY.out.versions)
    }
    if (diamond_enabled) {
        PREPARE_DIAMOND_QUERY(PREDICT_LONGEST_ATG_ORFS.out.orfs)
        diamond_query = PREPARE_DIAMOND_QUERY.out.query
        DIAMOND_BLASTP(
            diamond_query
                .filter { _meta, _backend, _scope, _query, status -> status.text.contains('\tVALID') }
                .map { meta, backend, scope, query, _status -> tuple(meta, backend, scope, query) },
            database
        )
        EMIT_EMPTY_DIAMOND_HITS(
            diamond_query.filter { _meta, _backend, _scope, _query, status -> status.text.contains('\tEMPTY') }
        )
        diamond_hits = DIAMOND_BLASTP.out.hits.mix(EMIT_EMPTY_DIAMOND_HITS.out.hits)
        WRITE_DATABASE_IDENTITY(database, 'diamond_database')
        REPORT_COMBINED_DIAMOND_MODELS(diamond_hits)
        comparison_inputs = validated_models
            .map { meta, backend, scope, bed12, _manifest, _validation -> tuple(meta.run_id ?: meta.cohort_id ?: meta.id, backend, scope, meta, bed12) }
            .join(REPORT_COMBINED_DIAMOND_MODELS.out.report.map { meta, backend, scope, report -> tuple(meta.run_id ?: meta.cohort_id ?: meta.id, backend, scope, report) }, by: [0, 1, 2])
            .map { _id, backend, scope, meta, bed12, report -> tuple(meta, backend, scope, bed12, report) }
        COMPARE_CANDIDATE_MODELS(comparison_inputs)
        diamond_reports = REPORT_COMBINED_DIAMOND_MODELS.out.report
        comparisons = COMPARE_CANDIDATE_MODELS.out.comparison
        diamond_versions = DIAMOND_BLASTP.out.versions.mix(WRITE_DATABASE_IDENTITY.out.versions).mix(REPORT_COMBINED_DIAMOND_MODELS.out.versions).mix(COMPARE_CANDIDATE_MODELS.out.versions)
    }
    emit:
    report = diamond_reports.mix(comparisons)
    orfs = orf_products
    versions = EXTRACT_COMBINED_TRANSCRIPTS.out.versions
        .mix(PREDICT_LONGEST_ATG_ORFS.out.versions)
        .mix(database_versions)
        .mix(diamond_versions)
        .mix(orf_versions)
}
