nextflow.enable.dsl = 2

include { EXTRACT_COMBINED_TRANSCRIPTS } from '../modules/extract_transcripts.nf'
include { PREDICT_LONGEST_ATG_ORFS } from '../modules/predict_orfs.nf'
include { PREP_DIAMOND_DB } from '../../qc/modules/diamond/prep_diamond_db.nf'
include { DIAMOND_BLASTP } from '../../qc/modules/diamond/blastp.nf'
include { REPORT_COMBINED_DIAMOND_MODELS } from '../modules/report_diamond_models.nf'

workflow RUN_DIAMOND_QC {
    take:
    combined_bed
    reference
    predictor
    reporter

    main:
    // combined_bed: tuple val(meta), path(combined_models.bed)
    if (!params.diamond_reference_db && !params.diamond_reference_proteins)
        error 'Provide --diamond_reference_db or --diamond_reference_proteins when --run_diamond_validation is enabled'
    EXTRACT_COMBINED_TRANSCRIPTS(combined_bed, reference)
    PREDICT_LONGEST_ATG_ORFS(EXTRACT_COMBINED_TRANSCRIPTS.out.transcripts, predictor)

    if (params.diamond_reference_db) {
        diamond_db = channel.value(file(params.diamond_reference_db, checkIfExists: true))
    } else {
        PREP_DIAMOND_DB(file(params.diamond_reference_proteins, checkIfExists: true))
        diamond_db = PREP_DIAMOND_DB.out.db
    }
    diamond_hits = DIAMOND_BLASTP(PREDICT_LONGEST_ATG_ORFS.out.peptides, diamond_db)
    report_input = PREDICT_LONGEST_ATG_ORFS.out.manifest
        .map { meta, manifest -> tuple(meta.id, meta, manifest) }
        .join(diamond_hits.hits.map { meta, hits -> tuple(meta.id, meta, hits) })
        .map { _id, meta, manifest, _hits_meta, hits -> tuple(meta, manifest, hits) }
    REPORT_COMBINED_DIAMOND_MODELS(report_input, diamond_db, reporter)

    version_ch = PREDICT_LONGEST_ATG_ORFS.out.versions
        .mix(DIAMOND_BLASTP.out.versions)
        .mix(REPORT_COMBINED_DIAMOND_MODELS.out.versions)
    if (!params.diamond_reference_db)
        version_ch = PREP_DIAMOND_DB.out.versions.mix(version_ch)

    emit:
    report = REPORT_COMBINED_DIAMOND_MODELS.out.report
    summary = REPORT_COMBINED_DIAMOND_MODELS.out.summary
    versions = version_ch
}
