// Collapse overlapping long-read transcript alignments into non-redundant
// gene/transcript models. Emulates the logic of HiveTranscriptCoalescer
// but operates entirely on flat files (BAM in, GFF3 out).

process COLLAPSE_TRANSCRIPTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::pysam=0.22 conda-forge::intervaltree=3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.1--py311h4e2aac4_0' :
        'biocontainers/pysam:0.22.1--py311h4e2aac4_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val   min_overlap      // reciprocal exon overlap fraction to merge transcripts
    val   max_intron_size  // filter alignments with introns larger than this

    output:
    tuple val(meta), path("*.collapsed.gff3"), emit: gff3
    tuple val(meta), path("*.stats.tsv"),      emit: stats
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    collapse_long_reads.py \\
        --bam         ${bam} \\
        --out         ${prefix}.collapsed.gff3 \\
        --stats       ${prefix}.stats.tsv \\
        --min-overlap ${min_overlap} \\
        --max-intron  ${max_intron_size} \\
        --sample-id   ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.collapsed.gff3 ${prefix}.stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.22.1
        python: 3.11
    END_VERSIONS
    """
}
