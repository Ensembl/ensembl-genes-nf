#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'
include { INVENTORY_LONG_READS } from './subworkflows/inventory_long_reads.nf'
include { PREPARE_LONG_READS } from './subworkflows/prepare_long_reads.nf'
include { ALIGN_LONG_READS } from './subworkflows/align_long_reads.nf'
include { RUN_ISOQUANT } from './subworkflows/run_isoquant.nf'
include { RUN_FLAIR } from './subworkflows/run_flair.nf'
include { RUN_BAMBU } from './subworkflows/run_bambu.nf'
include { COLLAPSE_LONG_READ_MODELS } from './subworkflows/collapse_long_read_models.nf'
include { MERGE_LONG_READ_MODELS } from './subworkflows/merge_long_read_models.nf'
include { VALIDATE_COMBINED_MODELS; VALIDATE_COMBINED_MODELS as VALIDATE_ACCESSION_MODELS } from './subworkflows/validate_combined_models.nf'
include { RUN_DIAMOND_QC } from './subworkflows/run_diamond_qc.nf'
include { RUN_DIAMOND_ANNOTATIONS } from './subworkflows/run_diamond_annotations.nf'
include { COLLECT_LONG_READ_SOFTWARE_VERSIONS } from './modules/collect_software_versions.nf'
include { COMPARE_CANDIDATE_MODELS; COMPARE_CANDIDATE_MODELS as COMPARE_ACCESSION_CANDIDATE_MODELS } from './modules/compare_candidate_models.nf'

def validate_runtime_options(auto_approve_safe, run_diamond_validation, annotation_only) {
    if (!(params.model_backend in ['tama', 'stringtie2', 'stringtie3', 'tmerge', 'isoquant', 'flair', 'bambu', 'all']))
        error '--model_backend must be tama, stringtie2, stringtie3, tmerge, isoquant, flair, bambu, or all'
    if (!(params.backend_merge_mode in ['native', 'legacy_common']))
        error "--backend_merge_mode must be native or legacy_common"
    if (!(params.backend_failure_policy in ['fail_fast', 'continue']))
        error "--backend_failure_policy must be fail_fast or continue"
    if (!(params.diamond_scope in ['cohort', 'accession', 'both']))
        error '--diamond_scope must be cohort, accession, or both'
    if (params.backend_merge_mode == 'native' && params.merge_tool != 'tama')
        error '--merge_tool is legacy-only; use --backend_merge_mode legacy_common when selecting tmerge'
    if (params.stringtie_merge_scope != 'run_then_cohort')
        error '--stringtie_merge_scope currently supports only run_then_cohort'
    if (!annotation_only && !params.manifest && !params.approved_manifest && !params.taxon_id)
        error 'Provide --taxon_id for discovery, --manifest for inventory, or --approved_manifest for production processing'
    if (annotation_only && (params.manifest || params.approved_manifest || params.taxon_id))
        error '--diamond_annotation_manifest cannot be combined with --manifest, --approved_manifest, or --taxon_id'
    if (annotation_only && !params.reference_fasta)
        error 'Provide --reference_fasta with --diamond_annotation_manifest'
    if (annotation_only && !run_diamond_validation)
        error '--run_diamond_validation must be enabled with --diamond_annotation_manifest'
    if (params.manifest && params.approved_manifest)
        error '--manifest and --approved_manifest cannot be combined'
    if (params.taxon_id && (params.manifest || params.approved_manifest))
        error '--taxon_id cannot be combined with --manifest or --approved_manifest'
    if (params.approved_manifest && !params.reference_fasta)
        error 'Provide --reference_fasta with the reference FASTA for production processing'
    if (auto_approve_safe && !params.manifest && !params.taxon_id)
        error 'Provide --manifest or --taxon_id with --auto_approve_safe'
    if (auto_approve_safe && params.approved_manifest)
        error '--auto_approve_safe cannot be combined with --approved_manifest'
    if ((auto_approve_safe || params.approved_manifest || params.taxon_id) && !params.reference_fasta)
        error 'Provide --reference_fasta for processing or automatic selection'
    if ((params.approved_manifest || auto_approve_safe || params.taxon_id) && !params.fastq_cache_dir)
        error 'Provide --fastq_cache_dir for production processing or automatic selection'
    if (run_diamond_validation && !params.diamond_reference_db && !params.diamond_reference_proteins)
        error 'Provide --diamond_reference_db or --diamond_reference_proteins when --run_diamond_validation is enabled'
    def isoquant_enabled = params.model_backend in ['isoquant', 'all']
    if (isoquant_enabled && !params.reference_fasta)
        error 'IsoQuant requires --reference_fasta'
    if (isoquant_enabled && params.isoquant_mode == 'annotation_free' && params.isoquant_genedb)
        error '--isoquant_genedb must be omitted in annotation_free mode'
    if (isoquant_enabled && params.isoquant_mode == 'reference_guided' && !params.isoquant_genedb)
        error '--isoquant_genedb is required in reference_guided mode'
    if (isoquant_enabled && !(params.isoquant_data_type in ['pacbio_ccs', 'nanopore', 'assembly']))
        error '--isoquant_data_type must be pacbio_ccs, nanopore, or assembly'
    if (isoquant_enabled && !(params.isoquant_scope in ['accession', 'cohort', 'both']))
        error '--isoquant_scope must be accession, cohort, or both'
    if (isoquant_enabled && !(params.isoquant_analysis instanceof List || params.isoquant_analysis))
        error '--isoquant_analysis must contain at least one analysis'
    if (isoquant_enabled && !(params.isoquant_analysis instanceof List ? params.isoquant_analysis.contains('transcript_discovery') : params.isoquant_analysis.toString().split().contains('transcript_discovery')))
        error '--isoquant_analysis must include transcript_discovery'
    if (isoquant_enabled && params.isoquant_args && params.isoquant_args =~ /--(reference|bam|genedb|output|data_type|analysis|threads|large_output|read_group)(\s|=|$)/)
        error '--isoquant_args cannot override pipeline-owned reference, BAM, genedb, output, data_type, analysis, threads, large_output, or read_group arguments'
    def flair_enabled = params.model_backend in ['flair', 'all']
    if (flair_enabled && params.flair_args && params.flair_args =~ /--(genomealignedbam|output|threads|genome|gtf|manifest)(\s|=|$)/)
        error '--flair_args cannot override pipeline-owned BAM, output, threads, genome, annotation, or manifest arguments'
    def bambu_enabled = params.model_backend in ['bambu', 'all']
    if (bambu_enabled && !(params.bambu_scope in ['accession', 'cohort', 'both']))
        error '--bambu_scope must be accession, cohort, or both'
}

def boolean_param(value) {
    value == true || value?.toString()?.toLowerCase() == 'true'
}

workflow {
    validateParameters()
    auto_approve_safe = boolean_param(params.auto_approve_safe)
    run_diamond_validation = boolean_param(params.run_diamond_validation)
    annotation_only = params.diamond_annotation_manifest != null
    validate_runtime_options(auto_approve_safe, run_diamond_validation, annotation_only)
    log.info(paramsSummaryLog(workflow))

    if (annotation_only) {
        RUN_DIAMOND_ANNOTATIONS(
            channel.value(file(params.diamond_annotation_manifest, checkIfExists: true)),
            file(params.reference_fasta, checkIfExists: true)
        )
        return
    }

    processing = params.approved_manifest || auto_approve_safe || params.taxon_id

    if (params.approved_manifest) {
        approved_source = channel.value(file(params.approved_manifest, checkIfExists: true))
    } else {
        INVENTORY_LONG_READS(
            params.manifest ? file(params.manifest, checkIfExists: true) : channel.empty(),
            params.taxon_id,
            params.metadata_json ? file(params.metadata_json, checkIfExists: true) : null,
            params.discovery_cache_dir ?: "${params.outdir}/discovery_cache",
            params.metadata_cache_dir ?: "${params.outdir}/metadata_cache",
            auto_approve_safe || params.taxon_id
        )
        approved_source = INVENTORY_LONG_READS.out.approved
    }

    inventory_versions = params.approved_manifest ? channel.empty() : INVENTORY_LONG_READS.out.versions

    if (!processing) {
        COLLECT_LONG_READ_SOFTWARE_VERSIONS(inventory_versions.collect())
        return
    }

    reference_fasta = file(params.reference_fasta, checkIfExists: true)
    cache = file(params.fastq_cache_dir).toAbsolutePath().toString()
    PREPARE_LONG_READS(
        approved_source,
        cache
    )
    ALIGN_LONG_READS(PREPARE_LONG_READS.out.reads, reference_fasta)
    isoquant_enabled = params.model_backend in ['isoquant', 'all']
    isoquant_bed = channel.empty()
    isoquant_accession_bed = channel.empty()
    isoquant_versions = channel.empty()
    isoquant_reports = channel.empty()
    if (isoquant_enabled) {
        RUN_ISOQUANT(ALIGN_LONG_READS.out.bam, reference_fasta)
        isoquant_bed = RUN_ISOQUANT.out.bed
        isoquant_accession_bed = RUN_ISOQUANT.out.accession_bed
        isoquant_versions = RUN_ISOQUANT.out.versions
        isoquant_reports = RUN_ISOQUANT.out.reports
    }
    flair_enabled = params.model_backend in ['flair', 'all']
    flair_bed = channel.empty()
    flair_accession_bed = channel.empty()
    flair_versions = channel.empty()
    flair_reports = channel.empty()
    if (flair_enabled) {
        flair_skips = (params.skip_flair_accessions ?: '').split(',').collect { value -> value.trim() }.findAll { value -> value }.toSet()
        flair_bams = ALIGN_LONG_READS.out.bam.filter { meta, _bam, _bai -> !flair_skips.contains(meta.id) }
        RUN_FLAIR(flair_bams, reference_fasta)
        flair_bed = RUN_FLAIR.out.bed
        flair_accession_bed = RUN_FLAIR.out.accession_bed
        flair_versions = RUN_FLAIR.out.versions
        flair_reports = RUN_FLAIR.out.reports
    }
    bambu_enabled = params.model_backend in ['bambu', 'all']
    bambu_bed = channel.empty()
    bambu_accession_bed = channel.empty()
    bambu_versions = channel.empty()
    if (bambu_enabled) {
        RUN_BAMBU(ALIGN_LONG_READS.out.bam, reference_fasta)
        bambu_bed = RUN_BAMBU.out.bed
        bambu_accession_bed = RUN_BAMBU.out.accession_bed
        bambu_versions = RUN_BAMBU.out.versions
    }

    collapse_versions = channel.empty()
    if (!(params.model_backend in ['isoquant', 'flair', 'bambu'])) {
        COLLAPSE_LONG_READ_MODELS(
            ALIGN_LONG_READS.out.bam,
            reference_fasta
        )
        collapse_versions = COLLAPSE_LONG_READ_MODELS.out.versions
    }
    merged_bed = channel.empty()
    merged_accession_bed = channel.empty()
    merge_versions = channel.empty()
    if (!(params.model_backend in ['isoquant', 'flair', 'bambu'])) {
        MERGE_LONG_READ_MODELS(COLLAPSE_LONG_READ_MODELS.out.native_models)
        merged_bed = MERGE_LONG_READ_MODELS.out.bed
        merged_accession_bed = MERGE_LONG_READ_MODELS.out.accession_bed
        merge_versions = MERGE_LONG_READ_MODELS.out.versions
    }
    cohort_candidates = merged_bed.mix(isoquant_bed).mix(flair_bed).mix(bambu_bed)
    accession_candidates = merged_accession_bed.mix(isoquant_accession_bed).mix(flair_accession_bed).mix(bambu_accession_bed)
    VALIDATE_COMBINED_MODELS(cohort_candidates)
    VALIDATE_ACCESSION_MODELS(accession_candidates)
    combined_bed = VALIDATE_COMBINED_MODELS.out.bed
    accession_bed = VALIDATE_ACCESSION_MODELS.out.bed
    comparison_input = combined_bed
        .map { _meta, _backend, bed -> bed }
        .collect()
        .map { beds -> tuple([id: params.cohort_id, scope: 'cohort'], beds, 'cohort') }
    COMPARE_CANDIDATE_MODELS(comparison_input)
    accession_comparison_input = accession_bed
        .map { _meta, _backend, bed -> bed }
        .collect()
        .map { beds -> tuple([id: 'accessions', scope: 'accession'], beds, 'accession') }
    COMPARE_ACCESSION_CANDIDATE_MODELS(accession_comparison_input)

    diamond_versions = channel.empty()
    if (run_diamond_validation) {
        cohort_diamond_inputs = combined_bed.map { meta, backend, bed ->
            tuple(meta + [scope: 'cohort'], backend, bed)
        }
        accession_diamond_inputs = accession_bed.map { meta, backend, bed ->
            tuple(meta + [scope: 'accession'], backend, bed)
        }
        diamond_inputs = params.diamond_scope == 'cohort' ? cohort_diamond_inputs :
            params.diamond_scope == 'accession' ? accession_diamond_inputs : cohort_diamond_inputs.mix(accession_diamond_inputs)
        RUN_DIAMOND_QC(diamond_inputs, reference_fasta)
        diamond_versions = RUN_DIAMOND_QC.out.versions
    }
    all_versions = channel.empty()
        .mix(inventory_versions)
        .mix(PREPARE_LONG_READS.out.versions)
        .mix(ALIGN_LONG_READS.out.versions)
        .mix(collapse_versions)
        .mix(merge_versions)
        .mix(isoquant_versions)
        .mix(flair_versions)
        .mix(bambu_versions)
        .mix(VALIDATE_COMBINED_MODELS.out.versions)
        .mix(diamond_versions)
    COLLECT_LONG_READ_SOFTWARE_VERSIONS(all_versions.collect())
}
