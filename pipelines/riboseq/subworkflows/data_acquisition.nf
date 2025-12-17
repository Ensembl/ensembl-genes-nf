/*
 * DATA ACQUISITION SUBWORKFLOW
 * Handles fetching and collapsing of ribosome profiling data
 */

include { LOCATE } from '../modules/locate.nf'
include { FASTQ_DL } from '../modules/fastq_dl.nf'
include { FASTQC } from '../modules/fastqc.nf'
include { FIND_ADAPTERS } from '../modules/find_adapters.nf'
include { EXTRACT_RPFS } from '../modules/extract_rpfs.nf'
include { FASTP } from '../modules/fastp.nf'
include { BOWTIE_RRNA_FILTER } from '../modules/bowtie_rrna_filter.nf'
include { COLLAPSE_FASTQ as COLLAPSE_FASTQ_INITIAL } from '../modules/collapse_fastq.nf'
include { COLLAPSE_FASTQ as COLLAPSE_FASTQ_FINAL } from '../modules/collapse_fastq.nf'

workflow DATA_ACQUISITION {
    take:
    sample_sheet      // path: CSV file with Run,study_accession columns
    adapter_list      // path: adapter list file for FASTQC

    main:
    // Validate required parameters
    if (params.fetch == null) {
        error "params.fetch must be defined as true or false"
    }

    // Parse sample sheet and create samples channel
    samples_ch = Channel
        .fromPath(sample_sheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id: row.Run,
                study_accession: row.study_accession ?: 'unknown',
            ]
            [ meta, row.Run ]
        }

    // Locate existing collapsed reads or runs that need processing
    LOCATE(samples_ch)

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
            collapsed_reads = Channel.empty()
        } else {
            // Process only samples that need processing
            needs_processing = LOCATE.out.needs_processing
            collapsed_reads = LOCATE.out.collapsed_reads
        }

        // Download FastQ files
        FASTQ_DL(needs_processing)

        // Run FastQC for quality control and adapter detection
        FASTQC(
            FASTQ_DL.out.fastq,
            adapter_list
        )

        // Branch based on architecture detection setting
        if (params.use_architecture_detection) {
            // Masterful RPF Extraction (Detects + Trims + Reports)
            EXTRACT_RPFS(
                FASTQ_DL.out.fastq,
                file(params.star_index)
            )

            // Re-collapse after extraction
            COLLAPSE_FASTQ_FINAL(EXTRACT_RPFS.out.rpfs)
            newly_collapsed_reads = COLLAPSE_FASTQ_FINAL.out.collapsed_fasta
        } else {
            // Traditional adapter finding approach
            FIND_ADAPTERS(
                FASTQ_DL.out.fastq,
                FASTQC.out.txt
            )

            // Trim with found adapters
            FASTP(
                FASTQ_DL.out.fastq
                    .join(FIND_ADAPTERS.out.adapter_report)
            )

            // Filter rRNA contamination (if rRNA index provided)
            if (params.rrna_index) {
                BOWTIE_RRNA_FILTER(
                    FASTP.out.trimmed_fastq,
                    file(params.rrna_index)
                )
                filtered_fastq = BOWTIE_RRNA_FILTER.out.filtered_fastq
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

        newly_collapsed_reads = Channel.empty()
        collapsed_reads = LOCATE.out.collapsed_reads
    }

    // Combine existing collapsed reads with newly processed ones
    all_collapsed_reads = collapsed_reads
        .mix(newly_collapsed_reads)
        .groupTuple()
        .map { meta, files -> tuple(meta, files[0]) }  // Take the first file if there are duplicates

    emit:
    samples = all_collapsed_reads
}
