/*
 * ORGANISM_SETUP SUBWORKFLOW
 * Prepare all reference files needed for the riboseq pipeline
 * Supports: Ensembl (gget), custom URLs, and SILVA rRNA
 */

include { GGET_DOWNLOAD } from '../modules/organism_setup/gget_download.nf'
include { URL_DOWNLOAD } from '../modules/organism_setup/url_download.nf'
include { BUILD_STAR_INDEX } from '../modules/organism_setup/build_star_index.nf'
include { MAKE_TRANSCRIPTOME } from '../modules/organism_setup/make_transcriptome.nf'
include { BUILD_BOWTIE_INDEX as BUILD_BOWTIE_TRANSCRIPTOME } from '../modules/organism_setup/build_bowtie_index.nf'
include { BUILD_BOWTIE_INDEX as BUILD_BOWTIE_RRNA } from '../modules/organism_setup/build_bowtie_index.nf'
include { EXTRACT_RRNA } from '../modules/organism_setup/extract_rrna.nf'
include { DOWNLOAD_SILVA } from '../modules/organism_setup/download_silva.nf'
include { GENERATE_CHROM_SIZES } from '../modules/organism_setup/generate_chrom_sizes.nf'
include { RIBOMETRIC_PREPARE } from '../modules/organism_setup/ribometric_prepare.nf'
include { GENERATE_CONFIG } from '../modules/organism_setup/generate_config.nf'

workflow ORGANISM_SETUP {
    take:
    organism          // val: organism name (e.g., "homo_sapiens")
    download_method   // val: 'gget' or 'url'
    version           // val: Ensembl version or 'custom'
    rrna_source       // val: 'gtf' or 'silva'
    fasta_url         // val: genome FASTA URL (only for url method)
    gtf_url           // val: GTF URL (only for url method)
    silva_url         // val: SILVA rRNA URL (only for silva source)
    gget_which        // val: list of file types for gget (default: ['dna', 'gtf'])

    main:
    // Initialize version tracking
    ch_versions = Channel.empty()

    // Normalize organism name: lowercase and replace spaces with underscores
    def organism_normalized = organism.toLowerCase().replace(' ', '_')

    //
    // Step 1: Download genome and annotation
    //
    if (download_method == 'gget') {
        // Download from Ensembl using gget
        GGET_DOWNLOAD(
            organism_normalized,
            version,
            gget_which ?: ['dna', 'gtf']
        )
        genome_fasta = GGET_DOWNLOAD.out.genome_fasta
        genome_gtf = GGET_DOWNLOAD.out.genome_gtf
        ch_versions = ch_versions.mix(GGET_DOWNLOAD.out.versions)
    } else if (download_method == 'url') {
        // Download from custom URLs
        URL_DOWNLOAD(
            fasta_url,
            gtf_url,
            organism_normalized,
            version ?: 'custom'
        )
        genome_fasta = URL_DOWNLOAD.out.genome_fasta
        genome_gtf = URL_DOWNLOAD.out.genome_gtf
        ch_versions = ch_versions.mix(URL_DOWNLOAD.out.versions)
    } else {
        error "Invalid download_method: ${download_method}. Choose 'gget' or 'url'."
    }

    //
    // Step 2: Build STAR index
    //
    BUILD_STAR_INDEX(
        genome_fasta,
        genome_gtf,
        organism_normalized,
        version
    )
    ch_versions = ch_versions.mix(BUILD_STAR_INDEX.out.versions)

    //
    // Step 3: Make transcriptome
    //
    MAKE_TRANSCRIPTOME(
        genome_gtf,
        genome_fasta,
        organism_normalized,
        version
    )
    ch_versions = ch_versions.mix(MAKE_TRANSCRIPTOME.out.versions)

    //
    // Step 4: Build Bowtie index for transcriptome
    //
    BUILD_BOWTIE_TRANSCRIPTOME(
        MAKE_TRANSCRIPTOME.out.transcripts,
        'transcriptome',
        organism_normalized,
        version
    )
    ch_versions = ch_versions.mix(BUILD_BOWTIE_TRANSCRIPTOME.out.versions)

    //
    // Step 5: Generate chromosome sizes
    //
    GENERATE_CHROM_SIZES(
        genome_fasta,
        organism_normalized,
        version
    )
    ch_versions = ch_versions.mix(GENERATE_CHROM_SIZES.out.versions)

    //
    // Step 6: Prepare RiboMetric annotation
    //
    RIBOMETRIC_PREPARE(
        genome_gtf,
        genome_fasta,
        organism_normalized,
        version
    )
    ch_versions = ch_versions.mix(RIBOMETRIC_PREPARE.out.versions)

    //
    // Step 8: Get rRNA sequences (from GTF or SILVA)
    //
    if (rrna_source == 'gtf') {
        EXTRACT_RRNA(
            genome_gtf,
            genome_fasta,
            organism_normalized,
            version
        )
        rrna_fasta = EXTRACT_RRNA.out.rrna_fasta
        ch_versions = ch_versions.mix(EXTRACT_RRNA.out.versions)
    } else if (rrna_source == 'silva') {
        DOWNLOAD_SILVA(
            silva_url,
            organism_normalized,
            version
        )
        rrna_fasta = DOWNLOAD_SILVA.out.rrna_fasta
        ch_versions = ch_versions.mix(DOWNLOAD_SILVA.out.versions)
    } else {
        error "Invalid rrna_source: ${rrna_source}. Choose 'gtf' or 'silva'."
    }

    //
    // Step 9: Build Bowtie index for rRNA
    //
    BUILD_BOWTIE_RRNA(
        rrna_fasta,
        'rRNA',
        organism_normalized,
        version
    )
    ch_versions = ch_versions.mix(BUILD_BOWTIE_RRNA.out.versions)

    //
    // Step 10: Generate config file for pipeline
    //
    GENERATE_CONFIG(
        BUILD_STAR_INDEX.out.index,
        BUILD_BOWTIE_TRANSCRIPTOME.out.index,
        BUILD_BOWTIE_RRNA.out.index,
        genome_gtf,
        genome_fasta,
        GENERATE_CHROM_SIZES.out.chrom_sizes,
        RIBOMETRIC_PREPARE.out.ribometric_tsv,
        MAKE_TRANSCRIPTOME.out.transcripts,
        organism_normalized,
        version
    )
    ch_versions = ch_versions.mix(GENERATE_CONFIG.out.versions)

    emit:
    star_index        = BUILD_STAR_INDEX.out.index              // path: STAR index directory
    bowtie_index      = BUILD_BOWTIE_TRANSCRIPTOME.out.index    // path: Bowtie1 transcriptome index
    rrna_index        = BUILD_BOWTIE_RRNA.out.index             // path: Bowtie1 rRNA index
    gtf               = genome_gtf                               // path: GTF annotation
    fasta             = genome_fasta                             // path: Genome FASTA
    chrom_sizes       = GENERATE_CHROM_SIZES.out.chrom_sizes    // path: Chromosome sizes file
    ribometric_anno   = RIBOMETRIC_PREPARE.out.ribometric_tsv   // path: RiboMetric annotation TSV
    transcriptome     = MAKE_TRANSCRIPTOME.out.transcripts      // path: Transcriptome FASTA
    config            = GENERATE_CONFIG.out.config              // path: Generated params.config
    versions          = ch_versions                             // channel: versions
}
