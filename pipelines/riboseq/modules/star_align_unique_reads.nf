/*
 * STAR alignment for unique reads FASTA
 * Single-pass alignment of deduplicated reads
 */

process STAR_ALIGN_UNIQUE_READS {
    tag "unique_reads"
    label 'process_high'
    cpus 16
    memory '48.GB'
    time '12.h'
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1

    conda "bioconda::star=2.7.11a bioconda::samtools=1.19"
    container "oras://community.wave.seqera.io/library/samtools_star:1b5dd3ca5b761fb8"

    publishDir "${params.outdir}/global", mode: 'copy'

    input:
    path fasta           // unique_reads.fasta
    path star_index      // STAR genome index

    output:
    path "unique_reads.bam", emit: bam
    path "unique_reads.bam.bai", emit: bai
    path "unique_reads_Log.final.out", emit: log
    path "versions.yml", emit: versions

    when:
    (task.ext.when == null || task.ext.when) && (params.matrix_align_unique_reads == null || params.matrix_align_unique_reads)

    script:
    def mismatches = params.mismatches ?: 2
    def max_multimappers = params.max_multimappers ?: 100
    def alignment_type = params.alignment_type ?: 'EndToEnd'
    """
    # STAR alignment with FASTA input
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fasta} \\
        --outFileNamePrefix unique_reads_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFilterMismatchNmax ${mismatches} \\
        --outFilterMultimapNmax ${max_multimappers} \\
        --alignEndsType ${alignment_type} \\
        --outSAMattributes NH HI AS NM MD \\
        --outFilterType BySJout \\
        --alignIntronMax 1 \\
        --alignSJoverhangMin 8 \\
        --limitBAMsortRAM ${task.memory.toBytes() - 1000000000}

    # Rename output
    mv unique_reads_Aligned.sortedByCoord.out.bam unique_reads.bam

    # Index
    samtools index unique_reads.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch unique_reads.bam
    touch unique_reads.bam.bai
    touch unique_reads_Log.final.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: unknown
        samtools: unknown
    END_VERSIONS
    """
}

process STAR_ALIGN_UNIQUE_READS_PARTITIONED {
    tag "${partition}"
    label 'process_high'
    cpus 16
    memory '48.GB'
    time '12.h'
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1

    conda "bioconda::star=2.7.11a bioconda::samtools=1.19"
    container "oras://community.wave.seqera.io/library/samtools_star:1b5dd3ca5b761fb8"

    publishDir "${params.outdir}/global_partitioned/${partition}", mode: 'copy'

    input:
    tuple val(partition), path(fasta)
    path star_index

    output:
    tuple val(partition), path("unique_reads.${partition}.bam"), emit: bam
    tuple val(partition), path("unique_reads.${partition}.bam.bai"), emit: bai
    tuple val(partition), path("unique_reads.${partition}_Log.final.out"), emit: log
    path "versions.yml", emit: versions

    when:
    (task.ext.when == null || task.ext.when) && (params.matrix_align_unique_reads == null || params.matrix_align_unique_reads)

    script:
    def mismatches = params.mismatches ?: 2
    def max_multimappers = params.max_multimappers ?: 100
    def alignment_type = params.alignment_type ?: 'EndToEnd'
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fasta} \\
        --outFileNamePrefix unique_reads.${partition}_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFilterMismatchNmax ${mismatches} \\
        --outFilterMultimapNmax ${max_multimappers} \\
        --alignEndsType ${alignment_type} \\
        --outSAMattributes NH HI AS NM MD \\
        --outFilterType BySJout \\
        --alignIntronMax 1 \\
        --alignSJoverhangMin 8 \\
        --limitBAMsortRAM ${task.memory.toBytes() - 1000000000}

    mv unique_reads.${partition}_Aligned.sortedByCoord.out.bam unique_reads.${partition}.bam
    samtools index unique_reads.${partition}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch unique_reads.${partition}.bam
    touch unique_reads.${partition}.bam.bai
    touch unique_reads.${partition}_Log.final.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: unknown
        samtools: unknown
    END_VERSIONS
    """
}
