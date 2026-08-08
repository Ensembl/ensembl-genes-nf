process RIBOWALTZ {
    tag "${meta.id}"

    conda "bioconda::bioconductor-ribowaltz=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:2.0--r43hdfd78af_0' :
        'biocontainers/ribowaltz:2.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.cds_coverage_psite.tsv.gz"), optional: true, emit: cds_coverage
    tuple val(meta), path("offset_plot"), optional: true, emit: offset_plots
    tuple val(meta), path("*.psite_offset.tsv.gz"), optional: true, emit: psite_offsets
    tuple val(meta), path("*.psite.tsv.gz"), optional: true, emit: psite_table
    tuple val(meta), path("*nt_coverage_psite.tsv.gz"), optional: true, emit: nt_coverage
    tuple val(meta), path("*.codon_coverage_rpf.tsv.gz"), optional: true, emit: codon_rpf
    tuple val(meta), path("*.codon_coverage_psite.tsv.gz"), optional: true, emit: codon_psite
    tuple val(meta), path("ribowaltz_qc/*.pdf"), optional: true, emit: qc_plots
    tuple val(meta), path("*.best_offset.txt"), optional: true, emit: best_offset
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def start_nts = params.ribowaltz_exclude_start ?: 0
    def stop_nts = params.ribowaltz_exclude_stop ?: 0
    def flanking = params.ribowaltz_flanking ?: 6
    def extremity = params.ribowaltz_extremity ?: 'auto'
    def cl = params.ribowaltz_confidence_level ?: 99
    def utr5l = params.ribowaltz_utr5_length ?: 25
    def cdsl = params.ribowaltz_cds_length ?: 40
    def utr3l = params.ribowaltz_utr3_length ?: 25
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(riboWaltz))
    suppressPackageStartupMessages(library(dplyr))

    # Create folders
    dir.create("offset_plot", showWarnings = FALSE)
    dir.create("ribowaltz_qc", showWarnings = FALSE)

    # Create annotation data table from GTF
    annotation_dt <- create_annotation("${gtf}")

    # Read BAM file
    bam_name <- "${prefix}"
    names(bam_name) <- sub(".bam\$", "", basename("${transcriptome_bam}"))

    reads_list <- bamtolist(bamfolder = ".", annotation = annotation_dt, name_samples = bam_name)

    # Filter reads - no filtering by default, use all lengths
    filtered_list <- reads_list

    # Check for reads overlapping start codon
    start_overlap <- filtered_list[[1]][end5 <= cds_start & end3 >= cds_start]
    if (nrow(start_overlap) == 0) {
        stop("No reads overlapping start codon. Cannot calculate P-site offsets.")
    }

    # Calculate P-site offsets
    # The psite() function with txt=TRUE automatically creates the best_offset.txt file
    psite_offset <- psite(filtered_list, flanking = ${flanking}, extremity = "${extremity}", start = TRUE,
                         txt = TRUE, plot = TRUE, plot_format = "pdf",
                         txt_file = "${prefix}.ribowaltz_best_offset.txt")

    # Export full offset table for this sample
    data.table::fwrite(psite_offset, "${prefix}.psite_offset.tsv.gz", sep = "\\t")

    # Create simplified offset file for bam_to_bed (only length and offset columns)
    # This matches the format expected by bam_to_bed.py: length<tab>offset
    # RiboWaltz uses different column names depending on extremity used
    # Try to find the correct offset column (auto selects best extremity)
    offset_col <- if ("corrected_offset_from_5" %in% names(psite_offset)) {
        "corrected_offset_from_5"
    } else if ("corrected_offset_from_3" %in% names(psite_offset)) {
        "corrected_offset_from_3"
    } else if ("offset_from_5" %in% names(psite_offset)) {
        "offset_from_5"
    } else if ("offset_from_3" %in% names(psite_offset)) {
        "offset_from_3"
    } else {
        stop("Could not find offset column in psite output")
    }

    offset_simple <- psite_offset[, c("length", offset_col), with = FALSE]
    colnames(offset_simple) <- c("length", "offset")
    write.table(offset_simple, file = "${prefix}.best_offset.txt",
                sep = "\\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    # Move offset plots to folder
    pdf_files <- list.files(pattern = "^offset.*\\\\.pdf\$")
    if (length(pdf_files) > 0) {
        file.rename(pdf_files, file.path("offset_plot", pdf_files))
    }

    # Update reads with P-site information
    filtered_psite_list <- psite_info(filtered_list, psite_offset, site = "psite",
                                     fasta_genome = TRUE, refseq_sep = " ",
                                     fastapath = "${fasta}",
                                     gtfpath = "${gtf}")

    # Write P-site table
    psite_table <- filtered_psite_list[[1]]
    psite_table <- mutate(psite_table, sample = "${prefix}")
    data.table::fwrite(psite_table, file = "${prefix}.psite.tsv.gz",
                      sep = "\\t")

    # Generate QC plots
    sample_name <- names(filtered_psite_list)[1]

    # Read length distribution
    length_dist <- rlength_distr(reads_list, sample = sample_name,
                                multisamples = "average", cl = ${cl},
                                colour = "grey70")
    ggplot2::ggsave("ribowaltz_qc/${prefix}_length_distribution.pdf",
                   length_dist[["plot"]], dpi = 400)

    # Meta-heatmap
    ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = sample_name,
                              cl = ${cl}, utr5l = ${utr5l}, cdsl = ${cdsl}, utr3l = ${utr3l})
    ggplot2::ggsave("ribowaltz_qc/${prefix}_ends_heatmap.pdf",
                   ends_heatmap[[paste0("plot_", sample_name)]],
                   dpi = 400, width = 12, height = 8)

    # P-site region distribution
    psite_region <- region_psite(filtered_psite_list, annotation = annotation_dt,
                                sample = sample_name)
    ggplot2::ggsave("ribowaltz_qc/${prefix}_psite_region.pdf",
                   psite_region[["plot"]], dpi = 400, width = 10)

    # Frame distribution
    min_length <- as.integer(min(psite_offset[,"length"]))
    max_length <- as.integer(max(psite_offset[,"length"]))

    frames_stratified <- frame_psite_length(filtered_psite_list, region = "all",
                                           sample = sample_name,
                                           length_range = min_length:max_length,
                                           annotation = annotation_dt)
    ggplot2::ggsave("ribowaltz_qc/${prefix}_frames_stratified.pdf",
                   frames_stratified[[paste0("plot_", sample_name)]],
                   dpi = 600, height = 9, width = 12)

    frames <- frame_psite(filtered_psite_list, region = "all",
                         length_range = min_length:max_length,
                         sample = sample_name, annotation = annotation_dt,
                         colour = "grey70")
    ggplot2::ggsave("ribowaltz_qc/${prefix}_frames.pdf",
                   frames[[paste0("plot_", sample_name)]],
                   dpi = 600, height = 9, width = 9)

    # Metaprofile
    metaprofile <- metaprofile_psite(filtered_psite_list, annotation_dt,
                                    sample = sample_name,
                                    utr5l = ${utr5l}, cdsl = ${cdsl}, utr3l = ${utr3l},
                                    colour = "black")
    ggplot2::ggsave("ribowaltz_qc/${prefix}_metaprofile_psite.pdf",
                   metaprofile[[paste0("plot_", sample_name)]],
                   dpi = 400, width = 12, height = 6)

    # Calculate coverage
    # Codon coverage
    rpf_coverage <- codon_coverage(filtered_psite_list, annotation = annotation_dt,
                                  sample = sample_name, psite = FALSE)
    data.table::fwrite(rpf_coverage, "${prefix}.codon_coverage_rpf.tsv.gz", sep = "\\t")

    psite_coverage <- codon_coverage(filtered_psite_list, annotation = annotation_dt,
                                    sample = sample_name, psite = TRUE)
    data.table::fwrite(psite_coverage, "${prefix}.codon_coverage_psite.tsv.gz", sep = "\\t")

    # CDS coverage
    cds_coverage <- cds_coverage(filtered_psite_list, annotation = annotation_dt)
    cols <- c("transcript", "length_cds", sample_name)
    cds_coverage_subset <- cds_coverage[, ..cols]
    data.table::fwrite(cds_coverage_subset, "${prefix}.cds_coverage_psite.tsv.gz", sep = "\\t")

    # CDS coverage with window
    cds_window_coverage <- cds_coverage(filtered_psite_list, annotation = annotation_dt,
                                       start_nts = ${start_nts}, stop_nts = ${stop_nts})
    cds_window_subset <- cds_window_coverage[, ..cols]
    data.table::fwrite(cds_window_subset,
                      "${prefix}.cds_plus${start_nts}nt_minus${stop_nts}nt_coverage_psite.tsv.gz",
                      sep = "\\t")

    # Write versions
    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    ribowaltz: ", packageVersion("riboWaltz"))
    ), "versions.yml")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p offset_plot ribowaltz_qc
    touch ${prefix}.cds_coverage_psite.tsv.gz
    touch ${prefix}.psite_offset.tsv.gz
    touch ${prefix}.psite.tsv.gz
    touch ${prefix}.codon_coverage_rpf.tsv.gz
    touch ${prefix}.codon_coverage_psite.tsv.gz
    touch ${prefix}.best_offset.txt
    touch offset_plot/${prefix}_offset_plot.pdf
    touch ribowaltz_qc/${prefix}_rlength_dist.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribowaltz: 2.0
    END_VERSIONS
    """
}
