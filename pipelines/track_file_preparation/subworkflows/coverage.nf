include { BAMCOVERAGE } from '../modules/bam_coverage.nf'

workflow PREPARE_COVERAGE_TRACKS {
    take:
    coverage_entities // tuple [meta, bam, bai, chrom_sizes]

    main:
    coverage = BAMCOVERAGE(coverage_entities)

    emit:
    results = coverage.results
    versions = coverage.versions
}
