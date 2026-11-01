#!/usr/bin/env nextflow

include { validateParameters } from 'plugin/nf-schema'
include { PLAN_AGGREGATES } from './subworkflows/plan_aggregates.nf'
include { AGGREGATE_COVERAGE_TRACKS } from './subworkflows/coverage.nf'
include { AGGREGATE_SPLICE_TRACKS } from './subworkflows/splice.nf'
include { WRITE_AGGREGATION_REPORTS } from './subworkflows/reporting.nf'

workflow {
    validateParameters()
    if (!params.track_manifest) error '--track_manifest is required'
    if (!params.aggregation_manifest) error '--aggregation_manifest is required'
    if (!params.chrom_sizes) error '--chrom_sizes is required'

    def track_manifest = file(params.track_manifest, checkIfExists: true)
    def aggregation_manifest = file(params.aggregation_manifest, checkIfExists: true)
    def manifest_base_dir = params.manifest_base_dir ?
        file(params.manifest_base_dir, checkIfExists: true).toString() :
        track_manifest.toFile().absoluteFile.parentFile.toString()
    def chrom_sizes = file(params.chrom_sizes, checkIfExists: true)

    input_manifests = channel.of(tuple(
        track_manifest,
        aggregation_manifest,
        manifest_base_dir,
        params.track_types,
        params.allow_opaque_normalization,
    ))
    plan = PLAN_AGGREGATES(
        input_manifests,
        file("${projectDir}/bin/validate_aggregation_inputs.py", checkIfExists: true),
    )

    coverage = AGGREGATE_COVERAGE_TRACKS(
        plan.plan,
        chrom_sizes,
        file("${projectDir}/bin/write_aggregate_record.py", checkIfExists: true),
    )
    splice = AGGREGATE_SPLICE_TRACKS(
        plan.plan,
        chrom_sizes,
        file("${projectDir}/resources/splice_junction.as", checkIfExists: true),
        file("${projectDir}/bin/aggregate_splice_junctions.py", checkIfExists: true),
        file("${projectDir}/bin/write_aggregate_record.py", checkIfExists: true),
    )

    results = coverage.results.mix(splice.results)
    versions = plan.versions.mix(coverage.versions).mix(splice.versions)
    WRITE_AGGREGATION_REPORTS(
        results,
        versions,
        plan.plan,
        file("${projectDir}/bin/write_aggregate_reports.py", checkIfExists: true),
    )
}
