// STAR_ALIGN
// Align paired-end (or single-end) RNA-seq reads with STAR.
// Produces sorted BAM + SJ.out.tab splice junction file.

process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::star=2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    publishDir "${params.outdir}/star_align/${meta.id}", mode: 'copy', pattern: "*.{bam,bai,tab}"

    input:
    tuple val(meta), path(reads)   // reads: [ R1.fq.gz ] or [ R1.fq.gz, R2.fq.gz ]
    path  index_dir               // STAR genome index directory

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"),    emit: bam
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam.bai"), emit: bai
    tuple val(meta), path("*.SJ.out.tab"),                        emit: sj
    tuple val(meta), path("*.Log.final.out"),                     emit: log
    path  "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: meta.id
    def args      = task.ext.args   ?: ''
    def reads_arg = reads instanceof List ? reads.join(' ') : reads
    def pe_flag   = (reads instanceof List && reads.size() > 1) ? '' : '--readFilesIn'
    """
    STAR \\
        --runThreadN     ${task.cpus} \\
        --genomeDir      ${index_dir} \\
        --readFilesIn    ${reads_arg} \\
        --readFilesCommand zcat \\
        --outSAMtype     BAM SortedByCoordinate \\
        --outSAMstrandField intronMotif \\
        --outSAMattributes NH HI AS NM \\
        --outFilterMultimapNmax 10 \\
        --outFilterMismatchNmax 10 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFileNamePrefix ${prefix}. \\
        --runRNGseed       0 \\
        ${args}

    samtools index ${prefix}.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam.bai
    touch ${prefix}.SJ.out.tab
    printf 'STAR   Mapping speed, Million of reads per hour\t0\\n' > ${prefix}.Log.final.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.11b
        samtools: 1.21
    END_VERSIONS
    """
}
