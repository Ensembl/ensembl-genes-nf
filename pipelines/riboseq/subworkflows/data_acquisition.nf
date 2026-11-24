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
            // All getRPF trimming completes before optional contaminant removal.
            EXTRACT_RPFS(
                FASTQ_DL.out.fastq,
                star_index
            )

            // The normal getRPF path keeps its existing collapsed output. When
            // contaminant removal is requested, use the fully trimmed FASTQ,
            // filter it, collapse the survivors, and only then apply the
            // downstream 20-40 nt gate.
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

            FILTER_RPF_LENGTHS(getrpf_collapsed_reads)
            newly_collapsed_reads = FILTER_RPF_LENGTHS.out.collapsed_fasta
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
        }
    } else {
        if (params.force_fetch) {
            log.warn "fetch is set to false but force_fetch is true. This combination is not valid. No data will be fetched or processed."
        } else {
            log.info "fetch is set to false. No data will be fetched or processed."
        }

        newly_collapsed_reads = channel.empty()
        collapsed_reads = LOCATE.out.collapsed_reads
    }

    // Combine existing collapsed reads with newly processed ones
    all_collapsed_reads = collapsed_reads
        .mix(newly_collapsed_reads)
        .groupTuple()
        .map { meta, files -> tuple(meta, files[0]) }  // Take the first file if there are duplicates

    emit:
    samples = all_collapsed_reads
    expected = expected_samples
}
