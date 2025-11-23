process FILTER_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::pysam=0.23.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.23.3--py39hdd5828d_0' :
        'biocontainers/pysam:0.23.3--py39hdd5828d_0' }"

    publishDir "${params.outdir}/filtered_bams", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.unique_no_junction.bam"), emit: unique_no_junction
    tuple val(meta), path("${prefix}.unique_with_junction.bam"), emit: unique_with_junction
    tuple val(meta), path("${prefix}.multi_no_junction.bam"), emit: multi_no_junction
    tuple val(meta), path("${prefix}.multi_with_junction.bam"), emit: multi_with_junction
    tuple val(meta), path("${prefix}*.bam"), emit: all_bams
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def max_mm = params.max_multimappers ?: 10
    def min_mapq = params.min_mapq ?: 0

    """
    python3 $projectDir/bin/filter_bam.py \\
        --bam ${bam} \\
        --prefix ${prefix} \\
        --max-multimappers ${max_mm} \\
        --min-mapq ${min_mapq} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unique_no_junction.bam
    touch ${prefix}.unique_with_junction.bam
    touch ${prefix}.multi_no_junction.bam
    touch ${prefix}.multi_with_junction.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
        pysam: 0.23.3
    END_VERSIONS
    """
}
