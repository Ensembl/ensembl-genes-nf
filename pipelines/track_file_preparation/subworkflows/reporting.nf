include { WRITE_TRACK_REPORTS } from '../modules/write_reports.nf'

workflow ASSEMBLE_TRACK_SETS {
    take:
    result_channels
    version_channels
    report_script

    main:
    results = result_channels.collect()
    versions = version_channels.collect()
    report = WRITE_TRACK_REPORTS(results, report_script, versions)

    emit:
    track_manifest = report.track_manifest
    entity_status = report.entity_status
    exceptions = report.exceptions
    versions = report.versions
}
