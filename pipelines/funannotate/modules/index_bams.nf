process INDEX_BAMS {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1@sha256:bf80e07e650becfd084db1abde0fe932b50f990a07fa56421ea647b552b5a406'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: indexed

    script:
    """
    for bam in staged/bam_*.bam; do
        bai="\${bam}.bai"
        if [[ ! -s "\$bai" ]]; then
            samtools index -@ ${task.cpus} "\$bam" "\$bai"
        fi
    done
    """

    stub:
    """
    mkdir -p staged
    for bam in staged/bam_*.bam; do
        touch "\${bam}.bai"
    done
    """
}
