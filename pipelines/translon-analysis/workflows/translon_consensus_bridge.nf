include { RENAME_BED } from '../../translon-consensus/modules/rename_bed.nf'
include { CREATE_SAMPLESHEET } from '../../translon-consensus/modules/create_samplesheet.nf'
include { REPORT_CONSENSUS } from '../../translon-consensus/modules/report_consensus.nf'
include { GENERATE_HTML_REPORT } from '../../translon-consensus/modules/generate_report.nf'

workflow TRANSLON_CONSENSUS_BRIDGE {
    take:
    bed12
    gtf

    main:
    renamed = RENAME_BED(bed12)
    grouped = renamed
        .map { meta, file -> tuple(meta.id, file) }
        .groupTuple()
        .map { id, files -> tuple([id: id], files.sort { file -> file.name }) }

    CREATE_SAMPLESHEET(grouped)
    REPORT_CONSENSUS(
        CREATE_SAMPLESHEET.out,
        params.ucsc_session_url,
        gtf
    )

    all_results = REPORT_CONSENSUS.out.results
        .map { _meta, files -> files }
        .flatten()
        .collect()
    GENERATE_HTML_REPORT(all_results, params.ucsc_session_url)

    emit:
    results = REPORT_CONSENSUS.out.results
    report = GENERATE_HTML_REPORT.out.html_report
}
