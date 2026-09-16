#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_LONG_READ_MANIFEST } from './modules/validate_manifest.nf'
include { DISCOVER_LONG_READ_DATA } from './modules/discover_long_read.nf'
include { INSPECT_LONG_READ_MANIFEST } from './modules/inspect_manifest.nf'
include { RESOLVE_LONG_READ_METADATA } from './modules/resolve_metadata.nf'
include { AUTO_APPROVE_SELECTED_LONG_READS } from './modules/auto_approve_selected.nf'
include { PREPARE_LONG_READS } from './subworkflows/prepare_long_reads.nf'
include { ALIGN_LONG_READS } from './subworkflows/align_long_reads.nf'
include { COLLAPSE_LONG_READ_MODELS } from './subworkflows/collapse_long_read_models.nf'
include { MERGE_LONG_READ_MODELS } from './subworkflows/merge_long_read_models.nf'
include { VALIDATE_COMBINED_MODELS } from './subworkflows/validate_combined_models.nf'
include { RUN_DIAMOND_QC } from './subworkflows/run_diamond_qc.nf'

workflow {
    if (!params.manifest && !params.approved_manifest && !params.taxon_id) error 'Provide --taxon_id for discovery, --manifest for inventory, or --approved_manifest for production processing'
    if (params.approved_manifest && !params.reference_fasta) error 'Provide --reference_fasta with the reference FASTA for production processing'
    if (params.auto_approve_safe && !params.manifest) error 'Provide --manifest with --auto_approve_safe'
    if (params.auto_approve_safe && params.approved_manifest) error '--auto_approve_safe cannot be combined with --approved_manifest'
    if (params.auto_approve_safe && !params.reference_fasta) error 'Provide --reference_fasta with --auto_approve_safe'
    if (!(params.shard_mode in ['none', 'contig'])) error "params.shard_mode must be 'none' or 'contig'"
    if (!(params.secondary_mode in ['no', 'yes'])) error "params.secondary_mode must be 'no' or 'yes'"
    if (!params.fastq_cache_dir && params.approved_manifest) error 'Provide --fastq_cache_dir for production processing'
    if (!params.fastq_cache_dir && params.auto_approve_safe) error 'Provide --fastq_cache_dir with --auto_approve_safe'
    auto_select = params.taxon_id && !params.manifest && !params.approved_manifest
    if (auto_select && !params.reference_fasta) error 'Provide --reference_fasta when selecting transcriptomic runs automatically'
    if (auto_select && !params.fastq_cache_dir) error 'Provide --fastq_cache_dir when selecting transcriptomic runs automatically'

    validator = file("${baseDir}/bin/read_input_classification.py", checkIfExists: true)
    if (auto_select) {
        DISCOVER_LONG_READ_DATA(
            params.taxon_id,
            file("${baseDir}/bin/discover_long_read_data.py", checkIfExists: true),
            file("${baseDir}/bin/read_input_classification.py", checkIfExists: true),
            file("${baseDir}/bin/resolve_metadata.py", checkIfExists: true),
            params.discovery_cache_dir ?: "${params.outdir}/discovery_cache"
        )
        selected_manifest = DISCOVER_LONG_READ_DATA.out.results.map { discovery_dir -> file("${discovery_dir}/proposed_long_read_manifest.tsv") }
        VALIDATE_LONG_READ_MANIFEST(selected_manifest, file("${baseDir}/bin/normalise_manifest.py", checkIfExists: true))
        if (params.metadata_json) {
            inventory_metadata = file(params.metadata_json, checkIfExists: true)
        } else {
            metadata_cache = params.metadata_cache_dir ?: "${params.outdir}/metadata_cache"
            RESOLVE_LONG_READ_METADATA(
                VALIDATE_LONG_READ_MANIFEST.out.manifest,
                file("${baseDir}/bin/resolve_metadata.py", checkIfExists: true),
                metadata_cache
            )
            inventory_metadata = RESOLVE_LONG_READ_METADATA.out.metadata
        }
        INSPECT_LONG_READ_MANIFEST(VALIDATE_LONG_READ_MANIFEST.out.manifest, inventory_metadata, validator)
        AUTO_APPROVE_SELECTED_LONG_READS(
            INSPECT_LONG_READ_MANIFEST.out.reports,
            file("${baseDir}/bin/auto_approve_selected.py", checkIfExists: true)
        )
        approved_source = AUTO_APPROVE_SELECTED_LONG_READS.out.manifest
    } else if (!params.approved_manifest) {
        VALIDATE_LONG_READ_MANIFEST(file(params.manifest, checkIfExists: true), file("${baseDir}/bin/normalise_manifest.py", checkIfExists: true))
        if (params.metadata_json) {
            inventory_metadata = file(params.metadata_json, checkIfExists: true)
        } else {
            metadata_cache = params.metadata_cache_dir ?: "${params.outdir}/metadata_cache"
            RESOLVE_LONG_READ_METADATA(
                VALIDATE_LONG_READ_MANIFEST.out.manifest,
                file("${baseDir}/bin/resolve_metadata.py", checkIfExists: true),
                metadata_cache
            )
            inventory_metadata = RESOLVE_LONG_READ_METADATA.out.metadata
        }
        INSPECT_LONG_READ_MANIFEST(VALIDATE_LONG_READ_MANIFEST.out.manifest, inventory_metadata, validator)
        if (params.auto_approve_safe) {
            AUTO_APPROVE_SELECTED_LONG_READS(
                INSPECT_LONG_READ_MANIFEST.out.reports,
                file("${baseDir}/bin/auto_approve_selected.py", checkIfExists: true)
            )
            approved_source = AUTO_APPROVE_SELECTED_LONG_READS.out.manifest
        } else {
            return
        }
    } else {
        approved_source = channel.value(file(params.approved_manifest, checkIfExists: true))
    }

    reference_fasta = file(params.reference_fasta, checkIfExists: true)
    cache = file(params.fastq_cache_dir).toAbsolutePath().toString()
    PREPARE_LONG_READS(approved_source, validator, cache)
    ALIGN_LONG_READS(PREPARE_LONG_READS.out.reads, reference_fasta)
    COLLAPSE_LONG_READ_MODELS(ALIGN_LONG_READS.out.bam, reference_fasta)
    merge_input = COLLAPSE_LONG_READ_MODELS.out.beds.map { beds -> tuple(params.cohort_id, beds.sort { left, right -> left.name <=> right.name }) }
    MERGE_LONG_READ_MODELS(merge_input)
    VALIDATE_COMBINED_MODELS(MERGE_LONG_READ_MODELS.out.bed, params.cohort_id)
    combined_bed = VALIDATE_COMBINED_MODELS.out.bed

    if (params.run_diamond_validation) {
        RUN_DIAMOND_QC(
            combined_bed,
            reference_fasta,
            file("${baseDir}/bin/longest_atg_orf.py", checkIfExists: true),
            file("${baseDir}/bin/report_diamond_models.py", checkIfExists: true)
        )
    }
}
