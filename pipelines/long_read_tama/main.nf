#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters ; paramsSummaryLog } from 'plugin/nf-schema'
include { INVENTORY_LONG_READS } from './subworkflows/inventory_long_reads.nf'
include { PREPARE_LONG_READS } from './subworkflows/prepare_long_reads.nf'
include { ALIGN_LONG_READS } from './subworkflows/align_long_reads.nf'
include { COLLAPSE_LONG_READ_MODELS } from './subworkflows/collapse_long_read_models.nf'
include { MERGE_LONG_READ_MODELS } from './subworkflows/merge_long_read_models.nf'
include { VALIDATE_COMBINED_MODELS } from './subworkflows/validate_combined_models.nf'
include { RUN_DIAMOND_QC } from './subworkflows/run_diamond_qc.nf'
include { COLLECT_LONG_READ_SOFTWARE_VERSIONS } from './modules/collect_software_versions.nf'

def validate_runtime_options(auto_approve_safe, run_diamond_validation) {
    if (!params.manifest && !params.approved_manifest && !params.taxon_id)
        error 'Provide --taxon_id for discovery, --manifest for inventory, or --approved_manifest for production processing'
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
}

def boolean_param(value) {
    value == true || value?.toString()?.toLowerCase() == 'true'
}

workflow {
    validateParameters()
    auto_approve_safe = boolean_param(params.auto_approve_safe)
    run_diamond_validation = boolean_param(params.run_diamond_validation)
    validate_runtime_options(auto_approve_safe, run_diamond_validation)
    log.info(paramsSummaryLog(workflow))

    validator = file("${baseDir}/bin/read_input_classification.py", checkIfExists: true)
    auditor = file("${baseDir}/bin/audit_fastq.py", checkIfExists: true)
    downloader = file("${baseDir}/bin/download_fastq.py", checkIfExists: true)
    acquirer = file("${baseDir}/bin/acquire_bam.py", checkIfExists: true)
    workload_inspector = file("${baseDir}/bin/inspect_bam_workload.py", checkIfExists: true)
    contig_splitter = file("${baseDir}/bin/split_bam_by_contig.py", checkIfExists: true)
    bed_validator = file("${baseDir}/bin/validate_tama_bed.py", checkIfExists: true)
    merge_filelist_builder = file("${baseDir}/bin/prepare_tama_merge_filelist.py", checkIfExists: true)
    processing = params.approved_manifest || auto_approve_safe || params.taxon_id

    if (params.approved_manifest) {
        approved_source = channel.value(file(params.approved_manifest, checkIfExists: true))
    } else {
        INVENTORY_LONG_READS(
            params.manifest ? file(params.manifest, checkIfExists: true) : channel.empty(),
            params.taxon_id,
            file("${baseDir}/bin/discover_long_read_data.py", checkIfExists: true),
            file("${baseDir}/bin/read_input_classification.py", checkIfExists: true),
            file("${baseDir}/bin/resolve_metadata.py", checkIfExists: true),
            file("${baseDir}/bin/extract_manifest_accessions.py", checkIfExists: true),
            file("${baseDir}/bin/normalise_manifest.py", checkIfExists: true),
            file("${baseDir}/bin/auto_approve_selected.py", checkIfExists: true),
            validator,
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
        validator,
        cache,
        auditor,
        downloader,
        acquirer
    )
    ALIGN_LONG_READS(PREPARE_LONG_READS.out.reads, reference_fasta)
    COLLAPSE_LONG_READ_MODELS(
        ALIGN_LONG_READS.out.bam,
        reference_fasta,
        workload_inspector,
        contig_splitter,
        bed_validator,
        merge_filelist_builder
    )
    merge_input = COLLAPSE_LONG_READ_MODELS.out.beds.map { beds -> tuple(params.cohort_id, beds.sort { left, right -> left.name <=> right.name }) }
    MERGE_LONG_READ_MODELS(merge_input, merge_filelist_builder)
    VALIDATE_COMBINED_MODELS(MERGE_LONG_READ_MODELS.out.bed, params.cohort_id, bed_validator)
    combined_bed = VALIDATE_COMBINED_MODELS.out.bed

    if (run_diamond_validation) {
        RUN_DIAMOND_QC(
            combined_bed,
            reference_fasta,
            file("${baseDir}/bin/longest_atg_orf.py", checkIfExists: true),
            file("${baseDir}/bin/report_diamond_models.py", checkIfExists: true)
        )
    }

    diamond_versions = run_diamond_validation ? RUN_DIAMOND_QC.out.versions : channel.empty()
    all_versions = channel.empty()
        .mix(inventory_versions)
        .mix(PREPARE_LONG_READS.out.versions)
        .mix(ALIGN_LONG_READS.out.versions)
        .mix(COLLAPSE_LONG_READ_MODELS.out.versions)
        .mix(MERGE_LONG_READ_MODELS.out.versions)
        .mix(VALIDATE_COMBINED_MODELS.out.versions)
        .mix(diamond_versions)
    COLLECT_LONG_READ_SOFTWARE_VERSIONS(all_versions.collect())
}
