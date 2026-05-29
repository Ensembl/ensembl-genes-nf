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
    - RiboMetric and RiboWaltz analysis (offsets, QC, profiles)
    - BEDgraph and BigWig generation
    - Track hub generation (optional)
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'

// Validate parameters against schema
validateParameters()

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

include { COLLECT_STAR_LOG } from './modules/collect_star_log.nf'
include { COLLECT_RIBOMETRIC } from './modules/collect_ribometric.nf'
include { COLLECT_GETRPF_CLEAN } from './modules/collect_getrpf_clean.nf'
include { QC_GATE } from './modules/qc_gate.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

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

        log.info "Organism setup complete. Config file generated at: ${params.outdir}/organism_setup/${params.organism}/${params.ensembl_version ?: 'custom'}/riboseq_params.config"
        return  // Exit after organism setup
    }

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

            log.info "Reference indices built. Proceeding with data processing..."
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
        if (!params.fasta) {
            error "Reference genome FASTA (--fasta) is required for RiboWaltz analysis."
        }
        if (!params.chrom_sizes_file) {
            error "Chromosome sizes file (--chrom_sizes_file) is required."
        }
    }

    // Reference files setup - choose source based on use_organism_setup flag
    star_index_ch = use_organism_setup ?
        ORGANISM_SETUP.out.star_index :
        Channel.fromPath(params.star_index, checkIfExists: true).first()

    gtf_ch = use_organism_setup ?
        ORGANISM_SETUP.out.gtf :
        Channel.fromPath(params.gtf, checkIfExists: true).first()

    fasta_ch = use_organism_setup ?
        ORGANISM_SETUP.out.fasta :
        Channel.fromPath(params.fasta, checkIfExists: true).first()

    chrom_sizes_ch = use_organism_setup ?
        ORGANISM_SETUP.out.chrom_sizes :
        Channel.fromPath(params.chrom_sizes_file, checkIfExists: true).first()

    ribometric_anno_ch = use_organism_setup ?
        ORGANISM_SETUP.out.ribometric_anno :
        (params.ribometric_annotation ?
            Channel.fromPath(params.ribometric_annotation, checkIfExists: true).first() :
            Channel.value(file('NO_FILE')))

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
        // SUBWORKFLOW: Analysis - RiboMetric and RiboWaltz
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

        // STAR alignment metrics
        COLLECT_STAR_LOG( run_id, ALIGNMENT.out.logs )

        // getRPF cleanliness
        COLLECT_GETRPF_CLEAN( run_id, QUALITY_CONTROL.out.reports, QUALITY_CONTROL.out.rpf_checks )

        // RiboMetric metrics + artifacts: join JSON, CSV, and offsets by sample
        ribometric_triplet = ANALYSIS.out.ribometric_json
            .join(ANALYSIS.out.ribometric_csv)
            .join(ANALYSIS.out.offsets)

        COLLECT_RIBOMETRIC(
            run_id,
            ribometric_triplet.map { meta, j, c, off -> [meta, j, c, off] }
        )

        // Gate using RiboMetric offsets; filter passing lengths
        def rules_path = params.qc_rules ?: "${projectDir}/resources/qc_rules.default.yaml"
        QC_GATE( run_id, ANALYSIS.out.offsets, file(rules_path) )

        // Use filtered offsets for downstream processing
        def offsets_for_post = QC_GATE.out.filtered_offsets

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
                .filter { meta, bw -> meta.bam_type == chosen_type }
                .map { meta, bw -> [ [ id: meta.id ], bw ] }
                .groupTuple(by: 0)
                .map { meta, files -> [ meta, files ] }

            // Gate by QC decision: selected_for_translon must be true
            def selected_ids = QC_GATE.out.qc_json
                .map { meta, qcjson ->
                    def js = new groovy.json.JsonSlurper().parse(qcjson.toFile())
                    js.selected_for_translon ? meta.id : null
                }
                .filter { it != null }

            // Key join on sample id
            def keyed_bw = bigwigs_per_sample.map { meta, files -> [ meta.id, [meta, files] ] }
            def keyed_sel = selected_ids.map { id -> [ id, true ] }
            allowed = keyed_bw.join(keyed_sel).map { id, pair, _ -> pair }

            TRANSLONSCORER(
                allowed,
                gtf_ch,
                fasta_ch
            )
        }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
