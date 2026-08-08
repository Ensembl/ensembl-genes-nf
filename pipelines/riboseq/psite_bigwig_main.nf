#!/usr/bin/env nextflow

/*
 * PSITE BIGWIG PIPELINE
 * Standalone pipeline to generate P-site bigwig files from paired genome/transcriptome BAMs
 *
 * Usage:
 *   nextflow run psite_bigwig_main.nf \
 *     --genome_bam_dir /path/to/genome_bams \
 *     --transcriptome_bam_dir /path/to/transcriptome_bams \
 *     --ribometric_annotation /path/to/annotation.tsv \
 *     --chrom_sizes /path/to/chrom.sizes \
 *     --outdir results
 */

nextflow.enable.dsl = 2

// Default parameters
params.genome_bam_dir = null
params.transcriptome_bam_dir = null
params.ribometric_annotation = null
params.chrom_sizes = null
params.outdir = 'results'
params.filter_bams = true

// RiboMetric parameters (with defaults)
params.ribometric_offset_method = 'changepoint'
params.ribometric_sample_size = 10000000

// Filter BAM parameters
params.max_multimappers = 10
params.min_mapq = 0

// Include subworkflow
include { PSITE_BIGWIG } from './subworkflows/psite_bigwig.nf'

workflow {
    // Validate required parameters inside the entry workflow for strict DSL2.
    if (!params.genome_bam_dir) { error "Please provide --genome_bam_dir" }
    if (!params.transcriptome_bam_dir) { error "Please provide --transcriptome_bam_dir" }
    if (!params.ribometric_annotation) { error "Please provide --ribometric_annotation" }
    if (!params.chrom_sizes) { error "Please provide --chrom_sizes" }

    // Load genome BAMs with their index files
    genome_bam_ch = channel
        .fromPath("${params.genome_bam_dir}/*.bam")
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai")
            def meta = [id: bam.baseName.replaceAll(/\.sorted$/, '')]
            tuple(meta, bam, bai)
        }

    // Load transcriptome BAMs with their index files
    transcriptome_bam_ch = channel
        .fromPath("${params.transcriptome_bam_dir}/*.bam")
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai")
            def meta = [id: bam.baseName.replaceAll(/\.sorted$/, '').replaceAll(/\.toTranscriptome$/, '')]
            tuple(meta, bam, bai)
        }

    // Load reference files as value channels
    ribometric_annotation_ch = channel.value(file(params.ribometric_annotation, checkIfExists: true))
    chrom_sizes_ch = channel.value(file(params.chrom_sizes, checkIfExists: true))

    // Run the pipeline
    PSITE_BIGWIG(
        genome_bam_ch,
        transcriptome_bam_ch,
        ribometric_annotation_ch,
        chrom_sizes_ch
    )
    log.info """
╔═══════════════════════════════════════════════════════════════╗
║                    P-SITE BIGWIG PIPELINE                     ║
╠═══════════════════════════════════════════════════════════════╣
║  Genome BAM dir       : ${params.genome_bam_dir}
║  Transcriptome BAM dir: ${params.transcriptome_bam_dir}
║  RiboMetric annotation: ${params.ribometric_annotation}
║  Chrom sizes          : ${params.chrom_sizes}
║  Output directory     : ${params.outdir}
║  Filter BAMs          : ${params.filter_bams}
║  Offset method        : ${params.ribometric_offset_method}
╚═══════════════════════════════════════════════════════════════╝
    """
}
