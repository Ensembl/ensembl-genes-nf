#!/usr/bin/env nextflow
/*
========================================================================================
    RIBOSEQ PIPELINE
========================================================================================
    Ribosome profiling data processing pipeline
    - Organism setup (optional): Build reference indices
    - Data acquisition and read collapsing
    - Quality control
    - STAR alignment (genome + transcriptome)
    - RiboMetric analysis (QC and offsets), with RiboWaltz offset fallback
    - BEDgraph and BigWig generation
    - Track hub generation (optional)
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Typed parameter declarations are required for correct CLI coercion under
// Nextflow 26 strict syntax. Configuration defaults remain in nextflow.config.
params {
    run_organism_setup: Boolean = false
    auto_build_indices: Boolean = false
    generate_trackhub: Boolean = false
    fetch: Boolean = true
    force_fetch: Boolean = false
    getrpf_max_reads: Integer = 5000
    getrpf_preserve_umi: Boolean = false
    min_read_length: Integer = 20
    max_read_length: Integer = 0
    run_rrna_filter: Boolean = false
    bowtie_rrna_mismatches: Integer = 2
    bowtie_rrna_k: Integer = 1
    ribodetector_chunk_size: Integer = 256
    ribodetector_len: Integer = 100
    mismatches: Integer = 3
    allow_introns: Boolean = true
    max_multimappers: Integer = 10
    min_mapq: Integer = 0
    trim_front: Integer = 0
    filter_by_length: Boolean = false
    rpf_length_min: Integer = 26
    rpf_length_max: Integer = 34
    ribometric_sample_size: Integer = 10000000
    ribometric_use_ribowaltz_offsets: Boolean = false
    run_choros: Boolean = false
    choros_num_genes: Integer = 250
    choros_min_coverage: Integer = 5
    choros_min_nonzero: Integer = 100
    ribowaltz_exclude_start: Integer = 0
    ribowaltz_exclude_stop: Integer = 0
    ribowaltz_flanking: Integer = 6
    ribowaltz_confidence_level: Integer = 99
    ribowaltz_utr5_length: Integer = 25
    ribowaltz_cds_length: Integer = 40
    ribowaltz_utr3_length: Integer = 25
    normalize_coverage: Boolean = false
    enable_unique_reads_tracking: Boolean = false
    run_matrix_mode: Boolean = false
    matrix_align_unique_reads: Boolean = true
    matrix_chunk_size: Integer = 10000
    matrix_sparse_shard_rows: Integer = 5000000
    matrix_sparse_read_bucket_size: Integer = 100000
    matrix_tsv_chunk_size: Integer = 500000
    matrix_metadata_shard_rows: Integer = 1000000
    matrix_use_partitioning: Boolean = false
    matrix_partition_prefix_length: Integer = 4
    matrix_partition_stride: Integer = 1000000000
    matrix_partition_sparse_read_bucket_size: Integer = 50000000
    matrix_partition_qc_enabled: Boolean = true
    matrix_partition_qc_min_total_records: Integer = 1
    matrix_partition_qc_min_total_counts: Integer = 1
    emit_both_offsets: Boolean = false
    translonscorer_plot_range: Integer = 30
    run_translonscorer: Boolean = false
}

include { validateParameters } from 'plugin/nf-schema'

/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

include { ORGANISM_SETUP } from './subworkflows/organism_setup.nf'
include { DATA_ACQUISITION } from './subworkflows/data_acquisition.nf'
include { QUALITY_CONTROL } from './subworkflows/quality_control.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { ANALYSIS } from './subworkflows/analysis.nf'
include { POST_PROCESSING } from './subworkflows/post_processing.nf'
include { TRACKHUB_GENERATION } from './subworkflows/trackhub_generation.nf'
include { UNIQUE_READS_MATRIX } from './subworkflows/unique_reads_matrix.nf'
include { TRANSLONSCORER } from './modules/translonscorer.nf'
include { CHOROS } from './modules/choros.nf'

include { COLLECT_QC_METRICS } from './modules/collect_qc_metrics.nf'
include { QC_GATE } from './modules/qc_gate.nf'
include { IMPORT_QC_DB } from './modules/import_qc_db.nf'
include { SAMPLE_COMPLETENESS } from './modules/sample_completeness.nf'
include { COLLATE_VERSIONS } from './modules/collate_versions.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Validate parameters inside the entry workflow for strict DSL2 syntax.
    validateParameters()

    // All modules publish their runtime-generated versions.yml files to the
    // shared versions topic. Collate them once after the workflow completes.
    def topic_versions = channel.topic('versions')
    def version_bundle = topic_versions
        .map { version_file -> "---\n${version_file.text}" }
        .collectFile(name: 'versions_bundle.yml')
    COLLATE_VERSIONS(version_bundle)

    //
    // MODE 1: Standalone organism setup only
    // Run with --run_organism_setup to only build reference files (no data processing)
    //
    if (params.run_organism_setup && !params.sample_sheet) {
        log.info "Running in ORGANISM SETUP mode (standalone)"

        // Validate required params for organism setup
        if (!params.organism) {
            error "Organism name (--organism) is required for organism setup"
        }
        if (!params.download_method) {
            error "Download method (--download_method) is required. Choose 'gget' or 'url'"
        }

        // Run organism setup
        ORGANISM_SETUP(
            params.organism,
            params.download_method,
            params.ensembl_version ?: 'custom',
            params.rrna_source ?: 'gtf',
            params.genome_fasta_url ?: '',
            params.genome_gtf_url ?: '',
            params.silva_url ?: '',
            params.gget_which ?: ['dna', 'gtf']
        )

        log.info "Organism setup complete. Reusable config: ${params.outdir}/pipeline_info/riboseq_params.config"
    } else {

    //
    // MODE 2: Full pipeline (with optional organism setup)
    //

    // Validate sample_sheet is provided for processing mode
        if (!params.sample_sheet) {
            error "Sample sheet (--sample_sheet) is required for data processing. Use --run_organism_setup without --sample_sheet to only build references."
        }

        //
        // Determine reference file sources
        // Priority: 1) Run organism setup if requested, 2) Use provided params
        //
        def use_organism_setup = false

        if (params.run_organism_setup || params.auto_build_indices) {
        // Check if we need to build indices
        def need_to_build = params.run_organism_setup ||
                           (params.auto_build_indices && (!params.star_index || !file(params.star_index).exists()))

        if (need_to_build) {
            use_organism_setup = true
            log.info "Building reference indices for ${params.organism}"

            // Validate required params
            if (!params.organism) {
                error "Organism name (--organism) is required when building indices"
            }

            ORGANISM_SETUP(
                params.organism,
                params.download_method ?: 'gget',
                params.ensembl_version ?: '115',
                params.rrna_source ?: 'gtf',
                params.genome_fasta_url ?: '',
                params.genome_gtf_url ?: '',
                params.silva_url ?: '',
                params.gget_which ?: ['dna', 'gtf']
            )

            log.info "Reference indices built. Reusable config: ${params.outdir}/pipeline_info/riboseq_params.config"
        }
        }

    // Validate params if not using organism setup
        if (!use_organism_setup) {
        if (!params.star_index) {
            error "STAR index (--star_index) is required. Use --run_organism_setup to build indices first."
        }
        if (!params.gtf) {
            error "GTF annotation (--gtf) is required."
        }
        if (!params.fasta && !params.ribometric_annotation) {
            error "Reference genome FASTA (--fasta) is required when --ribometric_annotation is not supplied, because the pipeline must fall back to RiboWaltz offsets."
        }
        if (!params.chrom_sizes_file) {
            error "Chromosome sizes file (--chrom_sizes_file) is required."
        }
        }

    // Reference files setup - choose source based on use_organism_setup flag
        star_index_ch = use_organism_setup ?
        ORGANISM_SETUP.out.star_index :
        channel.fromPath(params.star_index, checkIfExists: true).first()

        gtf_ch = use_organism_setup ?
        ORGANISM_SETUP.out.gtf :
        channel.fromPath(params.gtf, checkIfExists: true).first()

        fasta_ch = use_organism_setup ?
        ORGANISM_SETUP.out.fasta :
        (params.fasta ?
            channel.fromPath(params.fasta, checkIfExists: true).first() :
            channel.value(file('NO_FILE')))

        chrom_sizes_ch = use_organism_setup ?
        ORGANISM_SETUP.out.chrom_sizes :
        channel.fromPath(params.chrom_sizes_file, checkIfExists: true).first()

        ribometric_anno_ch = use_organism_setup ?
        ORGANISM_SETUP.out.ribometric_anno :
        (params.ribometric_annotation ?
            channel.fromPath(params.ribometric_annotation, checkIfExists: true).first() :
            channel.value(file('NO_FILE')))

        transcriptome_fasta_ch = use_organism_setup ?
        ORGANISM_SETUP.out.transcriptome :
        (params.transcriptome_fasta ?
            channel.fromPath(params.transcriptome_fasta, checkIfExists: true).first() :
            channel.value(file('NO_FILE')))

    if (params.run_choros) {
        if (!use_organism_setup && !params.ribometric_annotation) {
            error "--run_choros requires --ribometric_annotation (or organism setup)."
        }
        if (!use_organism_setup && !params.transcriptome_fasta) {
            error "--run_choros requires --transcriptome_fasta (or organism setup)."
        }
        if (!params.choros_container) {
            error "--run_choros requires a pinned --choros_container image. See containers/choros/Dockerfile."
        }
        if (params.ribometric_offset_target != 'a_site') {
            error "--run_choros requires --ribometric_offset_target a_site."
        }
    }

    //
    // SUBWORKFLOW: Data acquisition and read collapsing
    //
    DATA_ACQUISITION(
        params.sample_sheet,
        file(params.adapter_list),
        star_index_ch
    )

    //
    // SUBWORKFLOW: Quality control of collapsed reads
    //
    QUALITY_CONTROL(
        DATA_ACQUISITION.out.samples
    )

    //
    // MODE SELECTION: Matrix-based vs Per-sample analysis
    //
    // matrix_mode: Build global unique reads matrix, optionally align once, then downstream analysis
    // per_sample_mode: Traditional per-sample alignment with RiboMetric/RiboWaltz (default)
    //
    def use_matrix_mode = params.run_matrix_mode ?: false

    if (use_matrix_mode) {
        //
        // MATRIX MODE: Unique reads matrix pipeline
        // - Groups samples by study
        // - Builds study-level matrices
        // - Merges into global sparse counts by default, or dense Zarr/NPZ on request
        // - Optional single alignment of unique reads
        //
        log.info "Running in MATRIX MODE: Building global unique reads matrix"

        UNIQUE_READS_MATRIX(
            QUALITY_CONTROL.out.samples,
            star_index_ch
        )

        // Downstream analysis will use:
        // - UNIQUE_READS_MATRIX.out.global_counts / global_matrix
        // - UNIQUE_READS_MATRIX.out.unique_reads_bam (single BAM, if alignment enabled)
        // - UNIQUE_READS_MATRIX.out.global_metadata (read info)

        log.info "Matrix mode complete. Outputs in: ${params.outdir}/global/"

    } else {
        //
        // PER-SAMPLE MODE: Traditional per-sample alignment and analysis
        //
        log.info "Running in PER-SAMPLE MODE: Traditional alignment and analysis"

        //
        // SUBWORKFLOW: STAR alignment (genome + transcriptome)
        //
        ALIGNMENT(
            QUALITY_CONTROL.out.samples,
            star_index_ch,
            gtf_ch
        )

        //
        // SUBWORKFLOW: Analysis - prefer RiboMetric QC and offsets.
        // RiboWaltz is used only as an explicit or fallback offset source.
        //
        ANALYSIS(
            ALIGNMENT.out.transcriptome_bam,
            ribometric_anno_ch,
            gtf_ch,
            fasta_ch
        )

        //
        // SUBWORKFLOW: Post-processing (filter, BEDgraph, BigWig, unique reads index)
        //

        // Collect metrics into DuckDB and perform QC gating (RiboMetric-first)
        def run_id = workflow.runName

        // RiboMetric metrics + artifacts: join JSON, CSV, and offsets by sample
        ribometric_triplet = ANALYSIS.out.ribometric_json
            .join(ANALYSIS.out.ribometric_csv)
            .join(ANALYSIS.out.offsets)

        def rules_path = params.qc_rules ?: "${projectDir}/resources/qc_rules.default.yaml"

        def qc_metric_inputs = ALIGNMENT.out.logs
            .join(QUALITY_CONTROL.out.reports)
            .join(QUALITY_CONTROL.out.rpf_checks)
            .join(ribometric_triplet)
            .map { meta, star_log, getrpf_report, getrpf_checks, ribometric_json, ribometric_csv, offsets ->
                tuple(meta, star_log, getrpf_report, getrpf_checks, ribometric_json, ribometric_csv, offsets)
            }

        COLLECT_QC_METRICS(
            run_id,
            qc_metric_inputs
        )

        qc_gate_inputs = ANALYSIS.out.offsets
            .join(COLLECT_QC_METRICS.out.metrics)
            .map { meta, offsets, metrics ->
                tuple(meta, offsets, metrics)
            }

        QC_GATE(
            run_id,
            qc_gate_inputs,
            file(rules_path)
        )

        // Capture sample-level omissions caused by retries/ignored failures.
        // Structural validation errors still fail the workflow before this point.
        SAMPLE_COMPLETENESS(
            DATA_ACQUISITION.out.expected.collect(),
            DATA_ACQUISITION.out.samples.map { meta, _collapsed -> meta.id }.collect(),
            ALIGNMENT.out.transcriptome_bam.map { meta, _bam, _bai -> meta.id }.collect(),
            QC_GATE.out.qc_json.map { meta, _qc_json -> meta.id }.collect()
        )

        // Optional codon-level sequence-bias correction. Run once per sample
        // after the sample has passed the good QC gate; great is a stricter
        // subset and must not create a duplicate ChOROS task.
        if (params.run_choros) {
            choros_inputs = ALIGNMENT.out.transcriptome_bam
                .join(QC_GATE.out.good_offsets)
                .map { meta, bam, bai, offsets -> tuple(meta, bam, bai, offsets) }

            CHOROS(
                choros_inputs,
                ribometric_anno_ch,
                transcriptome_fasta_ch
            )
        }

        // Good and great QC tiers drive separate BEDgraph/BigWig tracks.
        // Failed samples still publish QC audit files, but do not enter
        // post-processing.
        def offsets_for_post = QC_GATE.out.good_offsets
                .map { meta, offsets -> tuple(meta + [track_tier: 'good'], offsets) }
            .mix(
                QC_GATE.out.great_offsets
                    .map { meta, offsets -> tuple(meta + [track_tier: 'great'], offsets) }
            )

        IMPORT_QC_DB(
            COLLECT_QC_METRICS.out.metrics.map { _meta, metrics -> metrics }.collect(),
            COLLECT_QC_METRICS.out.artifacts.map { _meta, artifacts -> artifacts }.collect(),
            QC_GATE.out.qc_rule_set.map { _meta, rule_set -> rule_set }.collect(),
            QC_GATE.out.qc_eval.map { _meta, qc_eval -> qc_eval }.collect(),
            QC_GATE.out.gate_selection.map { _meta, gate_selection -> gate_selection }.collect()
        )

        POST_PROCESSING(
            ALIGNMENT.out.genome_bam,
            offsets_for_post,
            chrom_sizes_ch,
            DATA_ACQUISITION.out.samples  // collapsed FASTA files
        )

        //
        // OPTIONAL: Track hub generation
        //
        if (params.generate_trackhub) {
            // Validate required params for trackhub
            if (!params.genome_assembly) {
                error "Genome assembly (--genome_assembly) is required for trackhub generation (e.g., 'hg38', 'mm39')"
            }

            TRACKHUB_GENERATION(
                POST_PROCESSING.out.bigwigs,
                POST_PROCESSING.out.merged_bigwigs,
                params.hub_name ?: 'riboseq_hub',
                params.genome_assembly,
                params.hub_email ?: 'noreply@example.com',
                params.sample_regex ?: '',
                params.annotation_regex ?: ''
            )
        }

        // OPTIONAL: TranslonScorer per-sample scoring (gated by QC)
        if (params.run_translonscorer) {
            // Select per-sample bigWigs for a chosen BAM type
            def chosen_type = params.translonscorer_bam_type ?: 'unique_no_junction'

            // Filter for the chosen BAM type and group multiple stranded bigwigs per sample
            bigwigs_per_sample = POST_PROCESSING.out.bigwigs
                .filter { meta, _bw -> meta.bam_type == chosen_type }
                .map { meta, bw -> tuple([ id: meta.id ], bw) }
                .groupTuple(by: 0)
                .map { meta, files -> tuple(meta, files) }

            // Gate by QC decision: selected_for_translon must be true
            def selected_ids = QC_GATE.out.translon_selected
                .map { meta, _selected -> meta.id }

            // Key join on sample id
            def keyed_bw = bigwigs_per_sample.map { meta, files -> tuple(meta.id, tuple(meta, files)) }
            def keyed_sel = selected_ids.map { id -> tuple(id, true) }
            allowed = keyed_bw.join(keyed_sel).map { _id, pair, _unused -> pair }

            TRANSLONSCORER(
                allowed,
                gtf_ch,
                fasta_ch
            )
        }
    }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
