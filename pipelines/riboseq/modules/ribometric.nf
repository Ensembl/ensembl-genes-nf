process RIBOMETRIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.10 conda-forge::biopython bioconda::pysam"
    container "ghcr.io/lapti-ucc/riboseqorg-nf-ribometric:latest"

    publishDir "${params.outdir}/RiboMetric", mode: 'copy'

    errorStrategy 'ignore'
    
    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
    path ribometric_annotation

    output:
    tuple val(meta), path("*RiboMetric.html"), emit: html
    tuple val(meta), path("*RiboMetric.json"), emit: json
    tuple val(meta), path("*RiboMetric.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RiboMetric run \\
        --bam ${transcriptome_bam} \\
        --annotation ${ribometric_annotation} \\
        --threads $task.cpus \\
        --html \\
        --json \\
        --csv \\
        --offset-calculation-method tripsviz \\
        -S 10000000 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: \$(RiboMetric --version 2>&1 | sed 's/RiboMetric version //g' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_RiboMetric.html
    touch ${prefix}_RiboMetric.json
    touch ${prefix}_RiboMetric.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: 1.0.0
    END_VERSIONS
    """
}
