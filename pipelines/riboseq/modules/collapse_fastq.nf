process COLLAPSE_FASTQ {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.10 pip::riboseq-dp-tools=0.1.10"
    container "ghcr.io/lapti-ucc/riboseqorg-nf-rdp-tools:latest"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*collapsed.fa"), emit: collapsed_fasta
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    RDP-Tools collapse $args $fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RDP-Tools: \$(RDP-Tools --version 2>&1 | head -n1 | sed 's/RDP-Tools //' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_collapsed.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RDP-Tools: unknown
    END_VERSIONS
    """
}
