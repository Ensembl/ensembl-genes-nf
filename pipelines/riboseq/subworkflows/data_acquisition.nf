/*
 * DATA ACQUISITION SUBWORKFLOW
 * Handles fetching and collapsing of ribosome profiling data
 *
 * Supports multiple adapter detection methods:
 * - getRPF (alignment-based extraction - RECOMMENDED)
 * - Traditional with multiple fastp modes:
 *   - fastp auto-detection
 *   - fastp with explicit adapter sequence
 *   - fastp with adapter FASTA from FASTQC top-hit
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
include { BOWTIE_RRNA_FILTER } from '../modules/bowtie_rrna_filter.nf'
include { RIBODETECTOR } from '../modules/ribodetector.nf'
include { COLLAPSE_FASTQ as COLLAPSE_FASTQ_FINAL } from '../modules/collapse_fastq.nf'

workflow DATA_ACQUISITION {
    take:
    sample_sheet      // path: CSV/TSV file with Run,study_accession columns
    adapter_list      // path: adapter list file for FASTQC
    star_index        // path: STAR index directory

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
            [ meta, row.Run ]
        }

    // Locate existing collapsed reads or runs that need processing
    LOCATE(samples_ch)

    // Preserve the complete expected run set so ignored sample-level failures
    // can be reported instead of disappearing from downstream joins.
    expected_samples = samples_ch.map { meta, id -> id }

    // Log runs that need processing
    LOCATE.out.needs_processing
        .map { meta, id, file -> id }
        .collectFile(name: "${params.outdir}/runs_needing_processing.txt", newLine: true)

    // Handle different combinations of fetch and force_fetch
    if (params.fetch) {
        if (params.force_fetch) {
            // Process all samples, including existing collapsed reads
            needs_processing = LOCATE.out.needs_processing
                .mix(
                    LOCATE.out.collapsed_reads.map { meta, collapsed_file ->
                        [meta, meta.id, collapsed_file]
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
            // Produces collapsed FASTA directly (no separate collapse step needed)
            EXTRACT_RPFS(
                FASTQ_DL.out.fastq,
                star_index
            )

            // Preserve all trimmed reads from getRPF, then gate the downstream RPF set.
            FILTER_RPF_LENGTHS(EXTRACT_RPFS.out.trimmed_collapsed)
            newly_collapsed_reads = FILTER_RPF_LENGTHS.out.collapsed_fasta
        } else {
            // Traditional fastp-based approach
            // Determine adapter source based on params:
            // 1. params.fastp_adapter_sequence - explicit sequence (handled in module)
            // 2. params.fastp_adapter_fasta - user-provided FASTA file
            // 3. FASTQC top-hit via FIND_ADAPTERS
            // 4. fastp auto-detection (no adapter file)

            if (params.fastp_adapter_fasta) {
                // User-provided adapter FASTA file
                adapter_file = channel.fromPath(params.fastp_adapter_fasta)
                fastq_with_adapter = FASTQ_DL.out.fastq
                    .combine(adapter_file)
            } else if (params.fastp_adapter_sequence) {
                // Explicit sequence - pass placeholder, module uses param
                fastq_with_adapter = FASTQ_DL.out.fastq
                    .map { meta, fastq -> [ meta, fastq, file('NO_ADAPTER_FILE') ] }
            } else if (params.adapter_detection_method == 'fastqc_tophit') {
                // Use FASTQC to find top-hit adapter
                FIND_ADAPTERS(
                    FASTQ_DL.out.fastq,
                    FASTQC.out.txt
                )
                fastq_with_adapter = FASTQ_DL.out.fastq
                    .join(FIND_ADAPTERS.out.adapter_report)
            } else {
                // fastp auto-detection (default for traditional)
                fastq_with_adapter = FASTQ_DL.out.fastq
                    .map { meta, fastq -> [ meta, fastq, file('NO_ADAPTER_FILE') ] }
            }

            // Unpack and run FASTP with appropriate inputs
            FASTP(
                fastq_with_adapter.map { meta, fastq, adapter -> [ meta, fastq ] },
                fastq_with_adapter.map { meta, fastq, adapter -> adapter }.first()
            )

            // Filter rRNA contamination (if explicitly enabled)
            if (params.run_rrna_filter) {
                def rrna_method = params.rrna_filter_method ?: 'bowtie'

                if (rrna_method == 'ribodetector') {
                    // ML-based rRNA detection (no reference required)
                    RIBODETECTOR(FASTP.out.trimmed_fastq)
                    filtered_fastq = RIBODETECTOR.out.filtered_fastq
                } else if (rrna_method == 'bowtie' && params.rrna_index) {
                    // Bowtie alignment-based filtering (requires index)
                    BOWTIE_RRNA_FILTER(
                        FASTP.out.trimmed_fastq,
                        file("${params.rrna_index}/*")
                    )
                    filtered_fastq = BOWTIE_RRNA_FILTER.out.filtered_fastq
                } else {
                    log.warn "rRNA filtering enabled but no valid method/index. Using unfiltered reads."
                    filtered_fastq = FASTP.out.trimmed_fastq
                }
            } else {
                filtered_fastq = FASTP.out.trimmed_fastq
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
