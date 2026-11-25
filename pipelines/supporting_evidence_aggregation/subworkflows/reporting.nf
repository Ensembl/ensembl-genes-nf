include { WRITE_SUPPORTING_EVIDENCE_REPORTS } from '../modules/write_reports.nf'

workflow WRITE_SUPPORTING_EVIDENCE_REPORTS_WORKFLOW {
    take:
    plan
    result_files
    writer

    main:
    report = WRITE_SUPPORTING_EVIDENCE_REPORTS(plan, result_files, writer)

    emit:
    manifest = report.manifest
    status = report.status
    metadata = report.metadata
    versions = report.versions
}
