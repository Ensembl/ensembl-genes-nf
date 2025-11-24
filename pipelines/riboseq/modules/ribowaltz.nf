process RIBOWALTZ {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bioconductor-ribowaltz=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:2.0--r43hdfd78af_0' :
        'biocontainers/ribowaltz:2.0--r43hdfd78af_0' }"

    publishDir "${params.outdir}/ribowaltz", mode: 'copy'

    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.cds_coverage_psite.tsv.gz"), optional: true, emit: cds_coverage
    tuple val(meta), path("offset_plot/*"), optional: true, emit: offset_plots
    tuple val(meta), path("*.psite_offset.tsv.gz"), optional: true, emit: psite_offsets
    tuple val(meta), path("*.psite.tsv.gz"), optional: true, emit: psite_table
    tuple val(meta), path("*nt_coverage_psite.tsv.gz"), optional: true, emit: nt_coverage
    tuple val(meta), path("*.codon_coverage_rpf.tsv.gz"), optional: true, emit: codon_rpf
    tuple val(meta), path("*.codon_coverage_psite.tsv.gz"), optional: true, emit: codon_psite
    tuple val(meta), path("ribowaltz_qc/*.pdf"), optional: true, emit: qc_plots
    tuple val(meta), path("*.best_offset.txt"), optional: true, emit: best_offset
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_start = params.ribowaltz_exclude_start ?: 0
    def exclude_stop = params.ribowaltz_exclude_stop ?: 0
    """
    #!/usr/bin/env Rscript

    library(riboWaltz)

    # Create annotation data table from GTF
    annotation_dt <- create_annotation(gtfpath = "${gtf}")

    # Read BAM file
    reads_list <- bamtolist(bamfolder = ".", annotation = annotation_dt)

    # Calculate P-site offsets
    psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")

    # Write offset table
    write.table(psite_offset, file = "${prefix}.psite_offset.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    system("gzip ${prefix}.psite_offset.tsv")

    # Write best offset summary
    best <- aggregate(percentage ~ sample, data = psite_offset, FUN = max)
    best <- merge(best, psite_offset)
    write.table(best[, c("sample", "length", "offset", "percentage")],
                file = "${prefix}.best_offset.txt",
                sep = "\\t", quote = FALSE, row.names = FALSE)

    # Create offset plots directory
    dir.create("offset_plot", showWarnings = FALSE)

    # Generate offset plots
    for (sample_name in names(reads_list)) {
        pdf(paste0("offset_plot/", sample_name, "_offset_plot.pdf"))
        print(psite_info(reads_list[[sample_name]], psite_offset))
        dev.off()
    }

    # Update reads with P-site information
    reads_psite_list <- psite_info(reads_list, psite_offset)

    # Write P-site table
    psite_table <- do.call(rbind, reads_psite_list)
    write.table(psite_table, file = "${prefix}.psite.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    system("gzip ${prefix}.psite.tsv")

    # Create QC directory
    dir.create("ribowaltz_qc", showWarnings = FALSE)

    # Generate QC plots
    for (sample_name in names(reads_psite_list)) {
        # Read length distribution
        pdf(paste0("ribowaltz_qc/", sample_name, "_rlength_dist.pdf"))
        print(rlength_distr(reads_list, sample_name))
        dev.off()

        # Read extremity heatmap
        pdf(paste0("ribowaltz_qc/", sample_name, "_rends_heat.pdf"))
        print(rends_heat(reads_list, annotation_dt, sample_name, cl = 85))
        dev.off()

        # Metaprofiles
        pdf(paste0("ribowaltz_qc/", sample_name, "_metaprofile.pdf"))
        print(metaprofile_psite(reads_psite_list, annotation_dt, sample_name))
        dev.off()

        # Frame distribution
        pdf(paste0("ribowaltz_qc/", sample_name, "_frame_psite.pdf"))
        print(frame_psite_length(reads_psite_list, sample_name))
        dev.off()
    }

    # Calculate CDS coverage
    cds_coverage <- cds_coverage(reads_psite_list, annotation_dt)
    write.table(cds_coverage, file = "${prefix}.cds_coverage_psite.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    system("gzip ${prefix}.cds_coverage_psite.tsv")

    # Calculate codon coverage
    codon_coverage_rpf <- codon_coverage(reads_list, annotation_dt, psite = FALSE)
    write.table(codon_coverage_rpf, file = "${prefix}.codon_coverage_rpf.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    system("gzip ${prefix}.codon_coverage_rpf.tsv")

    codon_coverage_psite <- codon_coverage(reads_psite_list, annotation_dt, psite = TRUE)
    write.table(codon_coverage_psite, file = "${prefix}.codon_coverage_psite.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    system("gzip ${prefix}.codon_coverage_psite.tsv")

    # Optional: calculate coverage excluding start/stop regions
    if (${exclude_start} > 0 || ${exclude_stop} > 0) {
        nt_coverage <- cds_coverage(reads_psite_list, annotation_dt,
                                   start_nts = ${exclude_start},
                                   stop_nts = ${exclude_stop})
        write.table(nt_coverage, file = "${prefix}.${exclude_start}nt_coverage_psite.tsv",
                   sep = "\\t", quote = FALSE, row.names = FALSE)
        system("gzip ${prefix}.${exclude_start}nt_coverage_psite.tsv")
    }

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
