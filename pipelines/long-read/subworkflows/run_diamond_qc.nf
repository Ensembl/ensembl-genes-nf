include { EXTRACT_COMBINED_TRANSCRIPTS; PREDICT_LONGEST_ATG_ORFS; PUBLISH_ORF_PRODUCTS; PREP_DIAMOND_DB; DIAMOND_BLASTP; REPORT_COMBINED_DIAMOND_MODELS; COMPARE_CANDIDATE_MODELS } from '../modules/finalization.nf'

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
        database_versions = PREP_DIAMOND_DB.out.versions
    }
    if (diamond_enabled) {
        DIAMOND_BLASTP(PREDICT_LONGEST_ATG_ORFS.out.orfs, database)
        REPORT_COMBINED_DIAMOND_MODELS(DIAMOND_BLASTP.out.hits)
        comparison_inputs = validated_models
            .map { meta, backend, scope, bed12, _manifest, _validation -> tuple(meta.run_id ?: meta.cohort_id ?: meta.id, backend, scope, meta, bed12) }
            .join(REPORT_COMBINED_DIAMOND_MODELS.out.report.map { meta, backend, scope, report -> tuple(meta.run_id ?: meta.cohort_id ?: meta.id, backend, scope, report) }, by: [0, 1, 2])
            .map { _id, backend, scope, meta, bed12, report -> tuple(meta, backend, scope, bed12, report) }
        COMPARE_CANDIDATE_MODELS(comparison_inputs)
        diamond_reports = REPORT_COMBINED_DIAMOND_MODELS.out.report
        comparisons = COMPARE_CANDIDATE_MODELS.out.comparison
        diamond_versions = DIAMOND_BLASTP.out.versions.mix(REPORT_COMBINED_DIAMOND_MODELS.out.versions).mix(COMPARE_CANDIDATE_MODELS.out.versions)
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
