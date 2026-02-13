#!/usr/bin/env nextflow

/*
 * PSITE BIGWIG PIPELINE
 * Standalone pipeline to generate P-site bigwig files from paired genome/transcriptome BAMs
 *
 * Usage:
 *   nextflow run psite_bigwig_main.nf \
 *     --genome_bam_dir /path/to/genome_bams \
 *     --transcriptome_bam_dir /path/to/transcriptome_bams \
 *     --gtf /path/to/annotation.gtf \
 *     --fasta /path/to/genome.fa \
 *     --chrom_sizes /path/to/chrom.sizes \
 *     --outdir results
 */

nextflow.enable.dsl = 2

// Default parameters
params.genome_bam_dir = null
params.transcriptome_bam_dir = null
params.gtf = null
params.fasta = null
params.chrom_sizes = null
params.outdir = 'results'
params.filter_bams = true

// RiboWaltz parameters (with defaults)
params.ribowaltz_exclude_start = 0
params.ribowaltz_exclude_stop = 0
params.ribowaltz_flanking = 6
params.ribowaltz_extremity = 'auto'
params.ribowaltz_confidence_level = 99
params.ribowaltz_utr5_length = 25
params.ribowaltz_cds_length = 40
params.ribowaltz_utr3_length = 25

// Filter BAM parameters
params.max_multimappers = 10
params.min_mapq = 0

// Validate required parameters
if (!params.genome_bam_dir) { error "Please provide --genome_bam_dir" }
if (!params.transcriptome_bam_dir) { error "Please provide --transcriptome_bam_dir" }
if (!params.gtf) { error "Please provide --gtf" }
if (!params.fasta) { error "Please provide --fasta" }
if (!params.chrom_sizes) { error "Please provide --chrom_sizes" }

// Include subworkflow
include { PSITE_BIGWIG } from './subworkflows/psite_bigwig.nf'

workflow {
    // Load genome BAMs with their index files
    genome_bam_ch = Channel
        .fromPath("${params.genome_bam_dir}/*.bam")
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai")
            def meta = [id: bam.baseName.replaceAll(/\.sorted$/, '')]
            [meta, bam, bai]
        }

    // Load transcriptome BAMs with their index files
    transcriptome_bam_ch = Channel
        .fromPath("${params.transcriptome_bam_dir}/*.bam")
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) bai = file("${bam.parent}/${bam.baseName}.bai")
            def meta = [id: bam.baseName.replaceAll(/\.sorted$/, '').replaceAll(/\.toTranscriptome$/, '')]
            [meta, bam, bai]
        }

    // Load reference files as value channels
    gtf_ch = Channel.value(file(params.gtf, checkIfExists: true))
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    chrom_sizes_ch = Channel.value(file(params.chrom_sizes, checkIfExists: true))

    // Run the pipeline
    PSITE_BIGWIG(
        genome_bam_ch,
        transcriptome_bam_ch,
        gtf_ch,
        fasta_ch,
        chrom_sizes_ch
    )
}

// Log pipeline info
log.info """
╔═══════════════════════════════════════════════════════════════╗
║                    P-SITE BIGWIG PIPELINE                     ║
╠═══════════════════════════════════════════════════════════════╣
║  Genome BAM dir      : ${params.genome_bam_dir}
║  Transcriptome BAM dir: ${params.transcriptome_bam_dir}
║  GTF                 : ${params.gtf}
║  FASTA               : ${params.fasta}
║  Chrom sizes         : ${params.chrom_sizes}
║  Output directory    : ${params.outdir}
║  Filter BAMs         : ${params.filter_bams}
╚═══════════════════════════════════════════════════════════════╝
"""
