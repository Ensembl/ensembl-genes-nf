process CONVERT_CRAMS {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1@sha256:bf80e07e650becfd084db1abde0fe932b50f990a07fa56421ea647b552b5a406'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: converted
    path 'versions.yml', emit: versions

    script:
    """
    genome=\$(find -L staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
    [[ -s "\$genome" ]]
    index=1
    for cram in staged/cram_*.cram; do
        [[ -s "\$cram" ]]
        bam="staged/bam_\$(printf '%02d' \$index).bam"
        samtools view -@ ${task.cpus} -T "\$genome" -b "\$cram" > "\$bam"
        [[ -s "\$bam" ]]
        index=\$((index + 1))
    done
    cat > versions.yml <<VERSIONS
"${task.process}":
    samtools: \$(samtools --version | head -n1)
VERSIONS
    """

    stub:
    """
    touch staged/bam_01.bam
    touch versions.yml
    """
}
