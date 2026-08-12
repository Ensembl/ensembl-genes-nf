include { ANNOTATION_INTEGRATE } from '../modules/annotation_integrate.nf'
include { SANITIZE_GFF } from '../modules/sanitize_gff.nf'

workflow ANNOTATION_REVIEW {
    take:
    sources
    integration_script
    sanitize_script

    main:
    SANITIZE_GFF(sources, sanitize_script)
    ANNOTATION_INTEGRATE(SANITIZE_GFF.out.sources, integration_script)

    emit:
    review = ANNOTATION_INTEGRATE.out.review
    cases = ANNOTATION_INTEGRATE.out.cases
    representative_cases = ANNOTATION_INTEGRATE.out.representative_cases
    havana_decisions = ANNOTATION_INTEGRATE.out.havana_decisions
    integrated = ANNOTATION_INTEGRATE.out.integrated
    session = ANNOTATION_INTEGRATE.out.session
    jbrowse_config = ANNOTATION_INTEGRATE.out.jbrowse_config
    validation_logs = ANNOTATION_INTEGRATE.out.validation_logs
    gffcompare = ANNOTATION_INTEGRATE.out.gffcompare
    readme = ANNOTATION_INTEGRATE.out.readme
    sanitization_reports = SANITIZE_GFF.out.reports
}
