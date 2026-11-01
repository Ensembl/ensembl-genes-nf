#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.enable.strict = true

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { PREPARE_LONG_READS } from './subworkflows/prepare_long_reads.nf'
include { ALIGN_LONG_READS } from './subworkflows/align_long_reads.nf'
include { PREPARE_CONTIG_SHARDS } from './subworkflows/prepare_contig_shards.nf'
include { COLLAPSE_TAMA } from './subworkflows/collapse_tama.nf'
include { COLLAPSE_STRINGTIE } from './subworkflows/collapse_stringtie.nf'
include { COLLAPSE_TMERGE } from './subworkflows/collapse_tmerge.nf'
include { RUN_ISOQUANT } from './subworkflows/run_isoquant.nf'
include { RUN_FLAIR } from './subworkflows/run_flair.nf'
include { RUN_BAMBU } from './subworkflows/run_bambu.nf'
include { MERGE_LONG_READ_MODELS } from './subworkflows/merge_long_read_models.nf'
include { NORMALIZE_VALIDATE_MODELS } from './subworkflows/normalize_validate_models.nf'
include { AUDIT_NATIVE_SHARDS } from './subworkflows/audit_native_shards.nf'
include { COMPARE_MODELS } from './subworkflows/compare_normalized_models.nf'
include { FINALIZE_LONG_READ_PIPELINE } from './subworkflows/evidence_finalize.nf'
include { COLLECT_LONG_READ_EVIDENCE_PRODUCTS } from './subworkflows/evidence_collection.nf'
include { RUN_DIAMOND_QC } from './subworkflows/run_diamond_qc.nf'

def selected_backends(value) {
    def supported = ['tama', 'stringtie2', 'stringtie3', 'tmerge', 'isoquant', 'flair', 'bambu']
    def requested = value instanceof List ? value : value.toString().split(',').collect { it.trim() }
    if (!requested || requested.size() != requested.toSet().size())
        error '--model_backends must be non-empty and contain no duplicates'
    if (requested == ['all']) return supported
    if (requested.contains('all') || requested.any { !supported.contains(it) })
        error "--model_backends values must be selected from ${supported.join(', ')} or the sole value all"
    requested
}

def validate_runtime_options() {
    if (!(params.mode in ['fastq', 'bam', 'approved_manifest']))
        error '--mode must be fastq, bam, or approved_manifest'
    if (!params.manifest && !params.approved_manifest)
        error 'an explicit --manifest or --approved_manifest is required'
    if (params.mode == 'approved_manifest' && !params.approved_manifest)
        error '--approved_manifest is required when --mode approved_manifest'
    if (params.mode != 'approved_manifest' && !params.manifest)
        error '--manifest is required for fastq and bam modes'
    if (!(params.ccs_mode in ['disabled', 'pbccs']))
        error '--ccs_mode must be disabled or pbccs'
    def stages = selected_preparation_stages(params.preparation_stages)
    if (stages.contains('bam_to_fastq') && params.mode == 'fastq')
        error '--preparation_stages bam_to_fastq requires bam or approved_manifest mode'
    if (stages.contains('bam_to_fastq') && params.ccs_mode == 'pbccs')
        error '--preparation_stages bam_to_fastq cannot be combined with ccs_mode pbccs'
    if (!params.reference_fasta && !params.reference_metadata)
        error '--reference_fasta or --reference_metadata is required'
    def backend_scopes = normalize_scope(params.backend_scopes)
    if (backend_scopes.isEmpty())
        error '--backend_scopes must contain run, cohort, or both'
    def backends = selected_backends(params.model_backends)
    if (backends.contains('isoquant') && !(params.isoquant_analysis instanceof List ? params.isoquant_analysis.contains('transcript_discovery') : params.isoquant_analysis.toString().split().contains('transcript_discovery')))
        error '--isoquant_analysis must include transcript_discovery when IsoQuant models are selected'
    if (backends.contains('isoquant') && params.isoquant_mode == 'reference_guided' && !params.isoquant_genedb)
        error '--isoquant_genedb is required for reference-guided IsoQuant'
    if (params.diamond_mode == 'enabled' && !params.diamond_database && !params.diamond_reference_proteins)
        error '--diamond_database or --diamond_reference_proteins is required when DIAMOND is enabled'
    reject_pipeline_owned_args('stringtie2_merge_args', params.stringtie2_merge_args, ['merge', 'output', 'gtf', 'reference'])
    reject_pipeline_owned_args('stringtie3_merge_args', params.stringtie3_merge_args, ['merge', 'output', 'gtf', 'reference'])
    reject_pipeline_owned_args('tmerge_args', params.tmerge_args, ['input', 'output', 'reference'])
    reject_pipeline_owned_args('isoquant_args', params.isoquant_args, ['bam', 'reference', 'genedb', 'output', 'data_type'])
    reject_pipeline_owned_args('isoquant_args', params.isoquant_args, ['large_output'])
    if (params.isoquant_large_output.contains('none') && params.isoquant_large_output.size() > 1)
        error '--isoquant_large_output none cannot be combined with other output types'
    reject_pipeline_owned_args('flair_args', params.flair_args, ['genomealignedbam', 'output', 'threads', 'genome', 'gtf', 'manifest', 'junction_bed', 'sample_name'])
    if (params.model_backends instanceof List && params.model_backends.contains('all') && params.model_backends.size() != 1)
        error 'all must be the only value in --model_backends'
    [backends: backends, backend_scopes: backend_scopes, preparation_stages: stages]
}

def normalize_scope(value) {
    def values = value instanceof List ? value : value.toString().split(',').collect { it.trim() }.findAll()
    if (values == ['both']) return ['run', 'cohort']
    if (values.any { !(it in ['run', 'cohort']) } || values.size() != values.toSet().size())
        error '--backend_scopes must contain only run and cohort without duplicates'
    values
}

def selected_preparation_stages(value) {
    def stages = value instanceof List ? value : value.toString().split(',').collect { it.trim() }.findAll()
    def supported = ['discover', 'resolve_metadata', 'acquire', 'bam_to_fastq']
    if (stages.size() != stages.toSet().size() || stages.any { !supported.contains(it) })
        error "--preparation_stages values must be selected from ${supported.join(', ')} without duplicates"
    stages
}

def reject_pipeline_owned_args(name, value, owned_flags) {
    def text = value?.toString()?.trim()
    def tokens = text ? text.split(/\s+/) : []
    if (text && owned_flags.any { flag -> tokens.any { token -> token == "--${flag}" || token.startsWith("--${flag}=") } })
        error "--${name} cannot override pipeline-owned ${owned_flags.collect { '--' + it }.join(', ')} arguments"
}

def approved_manifest_id(row) {
    row.id ?: row.run_accession ?: row.accession ?: row.sample_id
}

def normalized_read_representation(row, mode) {
    def raw = (row.read_representation ?: '').toString().trim().toLowerCase()
    def classification = (row.classification ?: '').toString().toUpperCase()
    def bam_input = row.bam || row.ccs?.toString()?.toLowerCase()?.endsWith('.bam') || classification.endsWith('_BAM')
    if (bam_input) {
        if (raw in ['raw_subread', 'raw_subreads'] || classification.contains('SUBREAD')) return 'raw_subreads'
        if (raw == 'ccs' || classification.contains('CCS')) return 'ccs'
        return 'bam'
    }
    mode == 'bam' && raw == 'bam' ? 'bam' : 'fastq'
}

def cached_approved_input(row, cache_dir) {
    if (!cache_dir || !row.selected_artifact_md5) return null
    def classification = row.classification?.toString()?.trim()
    def basename = row.selected_artifact_basename ?: row.input ?: row.fastq ?: row.bam ?: row.ccs
    if (!classification || !basename) return null
    def candidate = new File(cache_dir.toString(), "${classification}/${row.selected_artifact_md5}/${basename}")
    candidate.isFile() ? candidate.absolutePath : null
}

workflow {
    validateParameters()
    def runtime_options = validate_runtime_options()
    def backends = runtime_options.backends
    def backend_scopes = runtime_options.backend_scopes
    def preparation_stages = runtime_options.preparation_stages
    log.info(paramsSummaryLog(workflow))

    def manifest_path = params.mode == 'approved_manifest' ? params.approved_manifest : params.manifest
    manifest_file = file(manifest_path, checkIfExists: true)
    cache_dir = params.fastq_cache_dir ? file(params.fastq_cache_dir, checkIfExists: true).toAbsolutePath().toString() : null
    reference_path = params.reference_fasta
    if (!reference_path && params.reference_metadata) {
        def metadata_file = file(params.reference_metadata, checkIfExists: true)
        def metadata = new groovy.json.JsonSlurper().parse(metadata_file)
        reference_path = metadata.reference_fasta ?: metadata.fasta ?: metadata.path
        if (!reference_path) error '--reference_metadata must contain reference_fasta, fasta, or path'
        if (!(reference_path.toString().startsWith('/'))) reference_path = file("${metadata_file.parent}/${reference_path}").toString()
    }
    reference_file = file(reference_path, checkIfExists: true)
    if (params.reference_fai) file(params.reference_fai, checkIfExists: true)
    input_runs = channel.fromPath(manifest_file)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def row_id = approved_manifest_id(row)
            if (!row_id) error 'manifest row does not contain id, run_accession, accession, or sample_id'
            def source_value = row.source_uri ?: row.selected_artifact_uri ?: row.input ?: row.fastq ?: row.bam ?: row.ccs ?: row.path
            def cached_input = params.mode == 'approved_manifest' ? cached_approved_input(row, cache_dir) : null
            def input_name = cached_input ?: row.input ?: row.fastq ?: row.bam ?: row.ccs ?: row.path ?: source_value
            if (!input_name) error "manifest row ${row_id} does not contain an input path or source URI"
            if (!source_value) source_value = input_name
            def remote_source = source_value.startsWith('http://') || source_value.startsWith('https://') || source_value.startsWith('ftp://') || source_value.startsWith('s3://') || source_value.startsWith('gs://')
            def source_uri = remote_source ? source_value : (source_value.startsWith('/') ? source_value : file("${manifest_file.parent}/${source_value}").toString())
            def index_value = row.index_uri ?: row.bai ?: row.pbi
            def source_index_uri = index_value ? ((index_value.startsWith('/') || index_value.contains('://')) ? index_value : file("${manifest_file.parent}/${index_value}").toString()) : "${source_uri}.bai"
            def read_representation = normalized_read_representation(row, params.mode)
            def meta = [
                id: row_id,
                cohort_id: row.cohort_id ?: params.cohort_id,
                run_id: row.run_id ?: row_id,
                technology: row.technology ?: 'unknown',
                assembly_id: row.assembly_id ?: 'reference',
                read_set_id: row.read_set_id ?: row_id,
                sample_id: row.sample_id ?: row_id,
                library_id: row.library_id ?: row_id,
                read_representation: read_representation,
                source_read_representation: row.read_representation ?: '',
                classification: row.classification ?: '',
                selected_artifact_md5: row.selected_artifact_md5 ?: '',
                selected_artifact_basename: row.selected_artifact_basename ?: '',
                source_uri: source_uri,
                source_index_uri: source_index_uri
            ]
            if (preparation_stages.contains('acquire')) {
                tuple(meta, meta.source_uri)
            } else {
                def input_file = input_name.startsWith('/') ? file(input_name, checkIfExists: true) : file("${manifest_file.parent}/${input_name}", checkIfExists: true)
                if (meta.read_representation in ['bam', 'ccs', 'raw_subreads']) {
                if (meta.read_representation == 'raw_subreads' && params.ccs_mode != 'pbccs')
                    error "manifest row ${row_id} contains raw_subreads but --ccs_mode is not pbccs"
                def index_suffixes = meta.read_representation == 'raw_subreads' && params.ccs_mode == 'pbccs' ? ['.pbi', '.bai'] : ['.bai']
                def index_candidates = index_suffixes.collectMany { suffix -> [file("${input_file}${suffix}"), file(input_file.toString().replaceFirst(/\\.bam$/, suffix))] }
                def index_file = index_candidates.find { candidate -> candidate.exists() }
                if (!index_file) error "manifest row ${row_id} requires an indexed BAM sidecar"
                tuple(meta, input_file, index_file)
            } else {
                    tuple(meta, input_file)
                }
            }
        }

    PREPARE_LONG_READS(manifest_file, input_runs, params.mode, params.cohort_id, manifest_file.parent.toString(), params.ccs_mode, preparation_stages, params.discovery_query, cache_dir ?: '')
    ALIGN_LONG_READS(PREPARE_LONG_READS.out.reads, reference_file, params.mode, params.minimap2_preset, params.secondary_mode, params.alignment_qc)
    PREPARE_CONTIG_SHARDS(ALIGN_LONG_READS.out.alignment, params.shard_mode, params.shard_contig_reads, params.shard_very_large_reads)

    native_models = channel.empty()
    flair_combined_models = channel.empty()
    versions = PREPARE_LONG_READS.out.versions.mix(ALIGN_LONG_READS.out.versions).mix(PREPARE_CONTIG_SHARDS.out.versions)
    if (backends.contains('tama')) {
        COLLAPSE_TAMA(PREPARE_CONTIG_SHARDS.out.shards, reference_file, params.tama_cap_mode, params.tama_end_mode, params.tama_run_mode)
        native_models = native_models.mix(COLLAPSE_TAMA.out.native_models)
        versions = versions.mix(COLLAPSE_TAMA.out.versions)
    }
    if (backends.intersect(['stringtie2', 'stringtie3'])) {
        COLLAPSE_STRINGTIE(PREPARE_CONTIG_SHARDS.out.shards, reference_file, backends)
        native_models = native_models.mix(COLLAPSE_STRINGTIE.out.native_models)
        versions = versions.mix(COLLAPSE_STRINGTIE.out.versions)
    }
    if (backends.contains('tmerge')) {
        COLLAPSE_TMERGE(PREPARE_CONTIG_SHARDS.out.shards, params.tmerge_args)
        native_models = native_models.mix(COLLAPSE_TMERGE.out.native_models)
        versions = versions.mix(COLLAPSE_TMERGE.out.versions)
    }
    if (backends.contains('isoquant')) {
        RUN_ISOQUANT(ALIGN_LONG_READS.out.alignment, reference_file, params.isoquant_mode, params.isoquant_genedb ? file(params.isoquant_genedb, checkIfExists: true) : null, params.isoquant_data_type, params.isoquant_analysis, params.isoquant_high_memory, params.isoquant_args, params.isoquant_large_output)
        native_models = native_models.mix(RUN_ISOQUANT.out.native_models)
        versions = versions.mix(RUN_ISOQUANT.out.versions)
    }
    if (backends.contains('flair')) {
        RUN_FLAIR(PREPARE_CONTIG_SHARDS.out.shards, reference_file, params.flair_threads ?: 1, params.flair_args)
        native_models = native_models.mix(RUN_FLAIR.out.native_models)
        flair_combined_models = RUN_FLAIR.out.combined
        versions = versions.mix(RUN_FLAIR.out.versions)
    }
    if (backends.contains('bambu')) {
        RUN_BAMBU(ALIGN_LONG_READS.out.alignment, reference_file, params.bambu_track_reads, params.bambu_quantify)
        native_models = native_models.mix(RUN_BAMBU.out.native_models)
        versions = versions.mix(RUN_BAMBU.out.versions)
    }

    empty_native_inputs = PREPARE_CONTIG_SHARDS.out.summaries
        .filter { _meta, eligible_count -> eligible_count == 0 }
        .flatMap { meta, _eligible_count -> backends.collect { backend -> tuple(meta, [backend: backend, scope: 'run', run_id_or_cohort_id: meta.run_id, shard_id_or_merged: 'empty', native_format: 'gtf', model_count: 0, status: 'EMPTY']) } }
    native_models = native_models.mix(flair_combined_models)
    AUDIT_NATIVE_SHARDS(native_models, empty_native_inputs, PREPARE_CONTIG_SHARDS.out.plans, backends)
    audited_native_models = AUDIT_NATIVE_SHARDS.out.audited.map { meta, backend_meta, product, status, _audit -> tuple(meta, backend_meta, product, status) }
    versions = versions.mix(AUDIT_NATIVE_SHARDS.out.versions)

    scope_policy = backends.collectEntries { backend ->
        def selected = backend == 'isoquant' ? params.isoquant_scope : backend == 'bambu' ? params.bambu_scope : backend_scopes
        [(backend): normalize_scope(selected)]
    }
    merge_native_models = audited_native_models.filter { _meta, backend_meta, _product, _status -> backend_meta.backend != 'flair' || backend_meta.shard_id_or_merged in ['merged', 'empty'] }
    MERGE_LONG_READ_MODELS(merge_native_models, params.cohort_id, scope_policy, params.stringtie2_merge_args, params.stringtie3_merge_args, params.tmerge_args)
    versions = versions.mix(MERGE_LONG_READ_MODELS.out.versions)
    NORMALIZE_VALIDATE_MODELS(MERGE_LONG_READ_MODELS.out.merged, ALIGN_LONG_READS.out.reference_dictionary)
    versions = versions.mix(NORMALIZE_VALIDATE_MODELS.out.versions)
    comparison_records = NORMALIZE_VALIDATE_MODELS.out.validated.map { meta, backend, scope, bed12, _manifest, _validation -> [cohort_id: meta.cohort_id, backend: backend, scope: scope, model_path: bed12.name] }.collect()
    COMPARE_MODELS(comparison_records.map { groovy.json.JsonOutput.toJson(it) }, params.cohort_id)
    versions = versions.mix(COMPARE_MODELS.out.versions)
    evidence_status = channel.of('DISABLED')
    if (params.evidence_mode == 'enabled') {
        COLLECT_LONG_READ_EVIDENCE_PRODUCTS(NORMALIZE_VALIDATE_MODELS.out.validated, PREPARE_LONG_READS.out.normalized, ALIGN_LONG_READS.out.alignment_manifest.collect(), ALIGN_LONG_READS.out.reference_identity, params.cohort_id)
        evidence_status = COLLECT_LONG_READ_EVIDENCE_PRODUCTS.out.evidence.map { _evidence -> 'COMPLETE' }
        versions = versions.mix(COLLECT_LONG_READ_EVIDENCE_PRODUCTS.out.versions)
    }
    if (params.diamond_mode == 'enabled' || params.orf_mode == 'enabled') {
        diamond_db = params.diamond_database ? file(params.diamond_database, checkIfExists: true) : null
        diamond_proteins = params.diamond_reference_proteins ? file(params.diamond_reference_proteins, checkIfExists: true) : null
        RUN_DIAMOND_QC(NORMALIZE_VALIDATE_MODELS.out.validated, diamond_proteins, diamond_db, params.diamond_mode == 'enabled', params.orf_mode == 'enabled')
        versions = versions.mix(RUN_DIAMOND_QC.out.versions)
    }

    configuration = groovy.json.JsonOutput.toJson([
        pipeline_version: workflow.revision ?: 'uncommitted',
        nextflow_version: workflow.nextflow.version.toString(),
        mode: params.mode,
        manifest: manifest_file.toString(),
        reference_fasta: reference_file.toString(),
        reference_metadata: params.reference_metadata,
        reference_fai: params.reference_fai,
        cohort_id: params.cohort_id,
        fastq_cache_dir: params.fastq_cache_dir,
        model_backends: backends,
        requested_model_backends: params.model_backends,
        expanded_model_backends: backends,
        backend_scopes: backend_scopes,
        scope_policy: scope_policy,
        shard_mode: params.shard_mode,
        shard_contig_reads: params.shard_contig_reads,
        shard_very_large_reads: params.shard_very_large_reads,
        shard_max_concurrent: params.shard_max_concurrent,
        resource_retry_max: params.resource_retry_max,
        minimap2_preset: params.minimap2_preset,
        secondary_mode: params.secondary_mode,
        alignment_qc: params.alignment_qc,
        ccs_mode: params.ccs_mode,
        preparation_stages: preparation_stages,
        discovery_query: params.discovery_query,
        publish_native_intermediates: params.publish_native_intermediates,
        evidence_mode: params.evidence_mode,
        orf_mode: params.orf_mode,
        diamond_mode: params.diamond_mode,
        tama: [cap_mode: params.tama_cap_mode, end_mode: params.tama_end_mode, run_mode: params.tama_run_mode],
        stringtie2_merge_args: params.stringtie2_merge_args,
        stringtie3_merge_args: params.stringtie3_merge_args,
        tmerge_args: params.tmerge_args,
        isoquant: [mode: params.isoquant_mode, scope: params.isoquant_scope, data_type: params.isoquant_data_type, analysis: params.isoquant_analysis, high_memory: params.isoquant_high_memory, large_output: params.isoquant_large_output, args: params.isoquant_args],
        flair: [threads: params.flair_threads, args: params.flair_args],
        bambu: [scope: params.bambu_scope, track_reads: params.bambu_track_reads, quantify: params.bambu_quantify],
        diamond_database: params.diamond_database,
        diamond_reference_proteins: params.diamond_reference_proteins,
        tool_containers: params.tool_containers
    ])
    stage_status = groovy.json.JsonOutput.toJson([
        preparation: 'COMPLETE',
        alignment: 'COMPLETE',
        sharding: 'COMPLETE',
        reconstruction: backends.collectEntries { backend -> [(backend): 'COMPLETE'] },
        merge: 'COMPLETE',
        normalization: 'COMPLETE',
        validation: 'COMPLETE',
        orf: params.orf_mode == 'enabled' ? 'COMPLETE' : 'DISABLED',
        diamond: params.diamond_mode == 'enabled' ? 'COMPLETE' : 'DISABLED',
        evidence: params.evidence_mode == 'enabled' ? 'COMPLETE' : 'DISABLED',
        finalization: 'COMPLETE'
    ])
    normalized_pair = PREPARE_LONG_READS.out.normalized.map { normalized, checksum, input_checksums, status -> tuple(normalized, input_checksums, status) }
    validated_artifacts = NORMALIZE_VALIDATE_MODELS.out.validated.flatMap { _meta, _backend, _scope, bed12, manifest, validation -> [bed12, manifest, validation] }
    summary_artifacts = NORMALIZE_VALIDATE_MODELS.out.summary.map { _meta, _backend, _scope, summary -> summary }
    model_publications = NORMALIZE_VALIDATE_MODELS.out.validated.flatMap { meta, backend, scope, bed12, manifest, validation ->
        def identity = scope == 'run' ? meta.run_id : meta.cohort_id
        def root = scope == 'run' ? "per_run/${meta.run_id}" : "per_cohort/${meta.cohort_id}"
        [
            [source_path: bed12.name, published_path: "${root}/models/${backend}/normalized/${identity}.${backend}.bed12"],
            [source_path: manifest.name, published_path: "${root}/models/${backend}/normalized/${identity}.${backend}.model_manifest.json"],
            [source_path: validation.name, published_path: "${root}/models/${backend}/normalized/${identity}.${backend}.${scope}.model_validation.tsv"]
        ]
    }.mix(NORMALIZE_VALIDATE_MODELS.out.summary.map { meta, backend, scope, summary ->
        def identity = scope == 'run' ? meta.run_id : meta.cohort_id
        def root = scope == 'run' ? "per_run/${meta.run_id}" : "per_cohort/${meta.cohort_id}"
        [source_path: summary.name, published_path: "${root}/models/${backend}/normalized/${identity}.${backend}.${scope}.model_summary.tsv"]
    }).collect().map { records -> groovy.json.JsonOutput.toJson(records) }
    merge_provenance = MERGE_LONG_READ_MODELS.out.merged.map { meta, backend, scope, _gtf, filelist, filelist_checksum ->
        [backend: backend, scope: scope, run_id: meta.run_id, cohort_id: meta.cohort_id, inputs: filelist.text.readLines().findAll { it }, file_list_sha256: filelist_checksum.text.trim().split(/\s+/)[0]]
    }.collect().map { records -> groovy.json.JsonOutput.toJson(records) }
    FINALIZE_LONG_READ_PIPELINE(normalized_pair, PREPARE_CONTIG_SHARDS.out.plans.map { _meta, _bam, _bai, plan -> plan }.collect(), NORMALIZE_VALIDATE_MODELS.out.reports.mix(validated_artifacts).mix(summary_artifacts).collect(), reference_file, ALIGN_LONG_READS.out.reference_identity, versions.collect(), configuration, stage_status, evidence_status, model_publications, merge_provenance)
}
