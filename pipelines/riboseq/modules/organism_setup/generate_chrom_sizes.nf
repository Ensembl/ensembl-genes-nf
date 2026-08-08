/*
 * GENERATE_CHROM_SIZES
 * Generate chromosome sizes file from genome FASTA using samtools
 */

process GENERATE_CHROM_SIZES {
    tag "${organism}_${version}"
    label 'process_low'

    container "quay.io/biocontainers/samtools:0.1.19--2"

    input:
    path(genome_fasta)
    val(organism)
    val(version)

    output:
    path "*.chrom.sizes",   emit: chrom_sizes
    path "versions.yml",    emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = genome_fasta.baseName
    """
    # Generate .fai index file
    samtools faidx ${genome_fasta}

    # Extract chromosome sizes (columns 1 and 2 from .fai)
    cut -f1,2 ${genome_fasta}.fai > ${prefix}.chrom.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = genome_fasta.baseName
    """
    touch ${prefix}.chrom.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
