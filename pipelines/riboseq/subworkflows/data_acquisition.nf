/*
 * DATA ACQUISITION SUBWORKFLOW
 * Handles fetching and collapsing of ribosome profiling data
 *
 * Supports multiple adapter detection methods:
 * - getRPF (alignment-based extraction - RECOMMENDED)
 * - Traditional with selectable trimmers:
 *   - fastp
 *   - Trim Galore
 *   - cutadapt
 *
 * Supports multiple rRNA filtering methods:
 * - Bowtie (alignment-based, requires index)
 * - RiboDetector (ML-based, no reference required)
 */

include { LOCATE } from '../modules/locate.nf'
include { FASTQ_DL } from '../modules/fastq_dl.nf'
include { FASTQC } from '../modules/fastqc.nf'
include { FIND_ADAPTERS } from '../modules/find_adapters.nf'
include { EXTRACT_RPFS } from '../modules/extract_rpfs.nf'
include { FILTER_RPF_LENGTHS } from '../modules/filter_rpf_lengths.nf'
include { FASTP } from '../modules/fastp.nf'
include { TRIM_GALORE_RPF } from '../modules/trim_galore_rpf.nf'
include { CUTADAPT_RPF } from '../modules/cutadapt_rpf.nf'
include { BOWTIE_RRNA_FILTER } from '../modules/bowtie_rrna_filter.nf'
include { RIBODETECTOR } from '../modules/ribodetector.nf'
include { EXPAND_COLLAPSED_FASTQ } from '../modules/expand_collapsed_fastq.nf'
include { COLLAPSE_FASTQ as COLLAPSE_FASTQ_FINAL } from '../modules/collapse_fastq.nf'
include { COLLAPSE_FASTQ as COLLAPSE_FASTQ_GETRPF } from '../modules/collapse_fastq.nf'

workflow DATA_ACQUISITION {
    take:
    sample_sheet      // path: CSV/TSV file with Run,study_accession columns
    adapter_list      // path: adapter list file for FASTQC
    star_index        // path: STAR index directory
    rrna_index        // path: Bowtie1 rRNA index directory/files when enabled

    main:
    // Validate required parameters
    if (params.fetch == null) {
        error "params.fetch must be defined as true or false"
    }

    // Parse sample sheet and create samples channel
    // Sample sheet must have columns: Run, study_accession. Additional
    // metadata columns are allowed and ignored by this generic entrypoint.
    def sample_sheet_path = sample_sheet.toString()
    def sample_sheet_sep = params.sample_sheet_sep ?: (
        sample_sheet_path.toLowerCase().endsWith('.tsv') ? '\t' : ','
    )

    samples_ch = channel
        .fromPath(sample_sheet)
        .splitCsv(header: true, sep: sample_sheet_sep)
        .map { row ->
            def meta = [
                id: row.Run,
                study_id: row.study_accession ?: 'unknown',  // Used for matrix grouping
                study_accession: row.study_accession ?: 'unknown',  // Keep for backwards compat
            ]
            tuple(meta, row.Run)
        }

    // Normalize optional reported trims once, keyed by the raw accession so
    // every run receives its own recipe before library aggregation.
    if (params.getrpf_trim_recipes) {
        def recipe_path = params.getrpf_trim_recipes.toString()
        if (!file(params.getrpf_trim_recipes).exists()) error "getrpf_trim_recipes does not exist: ${recipe_path}"
        def recipe_sep = recipe_path.toLowerCase().endsWith('.tsv') ? '\t' : ','
        getrpf_recipe_map = channel.fromPath(params.getrpf_trim_recipes)
            .splitCsv(header: true, sep: recipe_sep)
            .collect()
            .map { rows ->
                def recipes = [:]
                def source_rows = [:]
                rows.eachWithIndex { row, index ->
                    def source = (row.source ?: "${recipe_path}:row${index + 2}").toString()
                    def run_id = row.run_id?.toString()?.trim()
                    if (!run_id) error "Malformed getRPF trim recipe at ${source}: run_id is blank"
                    def recipe = [run_id: run_id, source: source]
                    def adapter = row.adapter?.toString()?.trim()
                    recipe.adapter = adapter && adapter.toLowerCase() != 'passed' ? adapter : null
                    ['5p_barcode_length', '3p_barcode_length'].each { field ->
                        def raw = row[field]?.toString()?.trim()
                        if (!raw) {
                            recipe[field] = null
                        } else if (!(raw ==~ /[0-9]+/)) {
                            error "Malformed getRPF trim recipe at ${source}: ${field} must be a non-negative integer or blank"
                        } else {
                            recipe[field] = raw.toInteger()
                        }
                    }
                    def input_state = row.input_state?.toString()?.trim()
                    if (input_state && input_state != 'pretrimmed') {
                        error "Malformed getRPF trim recipe at ${source}: input_state must be blank or pretrimmed"
                    }
                    recipe.input_state = input_state ?: null
                    def comparable_recipe = recipe.findAll { key, _value -> key != 'source' }
                    def comparable_existing = recipes[run_id]?.findAll { key, _value -> key != 'source' }
                    if (recipes.containsKey(run_id) && comparable_existing != comparable_recipe) {
                        error "Conflicting getRPF trim recipes for ${run_id} at ${source_rows[run_id]} and ${source}"
                    }
                    recipes[run_id] = recipe
                    source_rows[run_id] = source
                }
                recipes
            }
    } else {
        getrpf_recipe_map = channel.value([:])
    }

    // Locate existing collapsed reads or runs that need processing
    LOCATE(samples_ch)

    // Preserve the complete expected run set so ignored sample-level failures
    // can be reported instead of disappearing from downstream joins.
    expected_samples = samples_ch.map { _meta, id -> id }

    // Log runs that need processing
    LOCATE.out.needs_processing
        .map { _meta, id, _file -> id }
        .collectFile(name: "${params.outdir}/runs_needing_processing.txt", newLine: true)

    // Handle different combinations of fetch and force_fetch
    if (params.fetch) {
        if (params.force_fetch) {
            // Process all samples, including existing collapsed reads
            needs_processing = LOCATE.out.needs_processing
                .mix(
                    LOCATE.out.collapsed_reads.map { meta, collapsed_file ->
                        tuple(meta, meta.id, collapsed_file)
                    }
                )
            collapsed_reads = channel.empty()
        } else {
            // Process only samples that need processing
            needs_processing = LOCATE.out.needs_processing
            collapsed_reads = LOCATE.out.collapsed_reads
        }

        // Download FastQ files
        FASTQ_DL(needs_processing)

        // Run FastQC for quality control (always useful for QC reports)
        FASTQC(
            FASTQ_DL.out.fastq,
            adapter_list
        )

        // Branch based on RPF extraction method
        // 'getrpf' = alignment-based extraction (recommended)
        // 'traditional' = fastp-based trimming + optional rRNA filter
        def use_getrpf = (params.rpf_extraction_method ?: 'getrpf') == 'getrpf'

        if (use_getrpf) {
            // getRPF Extraction (alignment-based - RECOMMENDED)
            // Uses STAR alignment to determine biological sequence boundaries
            // Reported trim recipes are paired with the raw accession before
            // getRPF runs, so every sample receives its own recipe.
            extraction_inputs = FASTQ_DL.out.fastq
                .combine(getrpf_recipe_map)
                .map { meta, fastq, recipes ->
                    tuple(meta, fastq, groovy.json.JsonOutput.toJson(recipes[meta.id]))
                }

            EXTRACT_RPFS(
                extraction_inputs,
                star_index
            )

            // getRPF has completed its trimming before this optional
            // contaminant-removal branch.  Filtering is performed on the
            // trimmed read stream, then survivors are collapsed again.
            if (params.run_rrna_filter) {
                def rrna_method = params.rrna_filter_method ?: 'bowtie'

                EXPAND_COLLAPSED_FASTQ(EXTRACT_RPFS.out.trimmed_collapsed)

                if (rrna_method == 'ribodetector') {
                    RIBODETECTOR(EXPAND_COLLAPSED_FASTQ.out.trimmed_fastq)
                    filtered_fastq = RIBODETECTOR.out.filtered_fastq
                } else if (rrna_method == 'bowtie') {
                    BOWTIE_RRNA_FILTER(
                        EXPAND_COLLAPSED_FASTQ.out.trimmed_fastq,
                        rrna_index
                    )
                    filtered_fastq = BOWTIE_RRNA_FILTER.out.filtered_fastq
                } else {
                    error "Unsupported rRNA filtering method: ${rrna_method}. Choose 'bowtie' or 'ribodetector'."
                }

                COLLAPSE_FASTQ_GETRPF(filtered_fastq)
                getrpf_collapsed_reads = COLLAPSE_FASTQ_GETRPF.out.collapsed_fasta
            } else {
                getrpf_collapsed_reads = EXTRACT_RPFS.out.trimmed_collapsed
            }

            if (params.getrpf_prepare_for_alignment) {
                // Read length is based on aligned query bases in preparation
                // mode, so defer the 20–40 nt policy until profile eligibility.
                newly_collapsed_reads = getrpf_collapsed_reads
                    .join(EXTRACT_RPFS.out.report)
                    .map { meta, reads, report ->
                        def report_json = new groovy.json.JsonSlurper().parse(report.toFile())
                        def local = (report_json.preparation ?: [:]).requires_local_alignment ?: false
                        def mode = local ? 'Local' : (params.alignment_type ?: 'EndToEnd')
                        tuple(meta + [alignment_mode: mode, preparation: true], reads)
                    }
            } else {
                FILTER_RPF_LENGTHS(getrpf_collapsed_reads)
                newly_collapsed_reads = FILTER_RPF_LENGTHS.out.collapsed_fasta
            }
            getrpf_hold_reads = params.align_getrpf_holds
                ? EXTRACT_RPFS.out.getrpf_hold_reads
                : channel.empty()
        } else {
            // Traditional adapter trimming.  All trimmers emit the same
            // FASTQ contract.  fastp and Trim Galore can auto-detect from the
            // reads directly; cutadapt receives a per-sample adapter FASTA
            // produced by FIND_ADAPTERS when no override is supplied.
            // Determine adapter source based on params:
            // 1. params.fastp_adapter_sequence - explicit sequence (handled in module)
            // 2. params.fastp_adapter_fasta - user-provided FASTA file
            // 3. FASTQC/adapter-list inference via FIND_ADAPTERS
            // 4. tool-specific default detection (no adapter file)

            def no_adapter_file = file("${projectDir}/resources/NO_ADAPTER_FILE")
            def traditional_trimmer = params.traditional_trimmer ?: 'fastp'
            def needs_pipeline_adapter_detection =
                params.adapter_detection_method == 'fastqc_tophit' ||
                (traditional_trimmer == 'cutadapt' &&
                    !params.fastp_adapter_sequence &&
                    !params.fastp_adapter_fasta)

            if (params.fastp_adapter_fasta) {
                // User-provided adapter FASTA file
                adapter_file = channel.fromPath(params.fastp_adapter_fasta, checkIfExists: true).first()
                fastq_with_adapter = FASTQ_DL.out.fastq
                    .combine(adapter_file)
            } else if (params.fastp_adapter_sequence) {
                // Explicit sequence is a parameter; pass the placeholder so
                // the tuple shape remains identical for every trimmer.
                fastq_with_adapter = FASTQ_DL.out.fastq.map { meta, fastq ->
                    tuple(meta, fastq, no_adapter_file)
                }
            } else if (needs_pipeline_adapter_detection) {
                // Use FASTQC to find top-hit adapter
                FIND_ADAPTERS(
                    FASTQ_DL.out.fastq,
                    FASTQC.out.txt,
                    file("${projectDir}/resources/adapter_list.tsv")
                )
                fastq_with_adapter = FASTQ_DL.out.fastq
                    .join(FIND_ADAPTERS.out.adapter_report)
            } else {
                // No external adapter file; fastp/Trim Galore use their
                // configured defaults.
                fastq_with_adapter = FASTQ_DL.out.fastq
                    .map { meta, fastq -> tuple(meta, fastq, no_adapter_file) }
            }

            if (traditional_trimmer == 'fastp') {
                FASTP(fastq_with_adapter)
                trimmed_fastq = FASTP.out.trimmed_fastq
            } else if (traditional_trimmer == 'trimgalore') {
                TRIM_GALORE_RPF(fastq_with_adapter)
                trimmed_fastq = TRIM_GALORE_RPF.out.trimmed_fastq
            } else if (traditional_trimmer == 'cutadapt') {
                CUTADAPT_RPF(fastq_with_adapter)
                trimmed_fastq = CUTADAPT_RPF.out.trimmed_fastq
            } else {
                error "Unsupported traditional_trimmer: ${traditional_trimmer}. Choose 'fastp', 'trimgalore', or 'cutadapt'."
            }

            // Filter rRNA contamination (if explicitly enabled)
            if (params.run_rrna_filter) {
                def rrna_method = params.rrna_filter_method ?: 'bowtie'

                if (rrna_method == 'ribodetector') {
                    // ML-based rRNA detection (no reference required)
                    RIBODETECTOR(trimmed_fastq)
                    filtered_fastq = RIBODETECTOR.out.filtered_fastq
                } else if (rrna_method == 'bowtie') {
                    // Bowtie alignment-based filtering (requires the wired index)
                    BOWTIE_RRNA_FILTER(
                        trimmed_fastq,
                        rrna_index
                    )
                    filtered_fastq = BOWTIE_RRNA_FILTER.out.filtered_fastq
                } else {
                    log.warn "rRNA filtering enabled but no valid method/index. Using unfiltered reads."
                    filtered_fastq = trimmed_fastq
                }
            } else {
                filtered_fastq = trimmed_fastq
            }

            // Collapse after trimming and filtering
            COLLAPSE_FASTQ_FINAL(filtered_fastq)
            newly_collapsed_reads = COLLAPSE_FASTQ_FINAL.out.collapsed_fasta
            getrpf_hold_reads = channel.empty()
        }
    } else {
        if (params.force_fetch) {
            log.warn "fetch is set to false but force_fetch is true. This combination is not valid. No data will be fetched or processed."
        } else {
            log.info "fetch is set to false. No data will be fetched or processed."
        }

        newly_collapsed_reads = channel.empty()
        collapsed_reads = LOCATE.out.collapsed_reads
        getrpf_hold_reads = channel.empty()
    }

    // Combine generations without allowing a reused output directory to pick
    // an arbitrary file. A duplicate sample is valid only when the lineage
    // keys agree; otherwise fail closed before downstream alignment.
    all_collapsed_reads = collapsed_reads
        .mix(newly_collapsed_reads)
        .map { meta, collapsed_file -> tuple(meta.id, meta, collapsed_file) }
        .groupTuple()
        .map { id, metas, files ->
            def lineage = metas.collect { meta ->
                meta.lineage_id ?: meta.run_id ?: meta.id
            }.unique()
            if (lineage.size() != 1) {
                error "Duplicate or conflicting collapsed reads for ${id}; refusing arbitrary file selection (lineage=${lineage})"
            }
            // The channel order is existing generation, then newly processed
            // generation. Matching lineage permits an explicit newest-wins
            // choice; conflicting lineage still fails closed.
            tuple(metas[-1], files[-1])
        }

    emit:
    samples = all_collapsed_reads
    expected = expected_samples
    getrpf_hold_reads = getrpf_hold_reads
}
