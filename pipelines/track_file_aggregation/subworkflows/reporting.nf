include { WRITE_AGGREGATE_REPORTS } from '../modules/write_aggregate_reports.nf'

workflow WRITE_AGGREGATION_REPORTS {
    take:
    results
    versions
    plan
    report_script

    main:
    result_files = results.map { _meta, _track, result, _provenance -> result }.collect()
    report = WRITE_AGGREGATE_REPORTS(result_files, versions.collect(), plan, report_script)

    emit:
    track_manifest = report.track_manifest
    aggregate_status = report.aggregate_status
    exceptions = report.exceptions
    versions = report.versions
}
