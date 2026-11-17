process CONVERT_CRAMS_TO_BAMS {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1@sha256:bf80e07e650becfd084db1abde0fe932b50f990a07fa56421ea647b552b5a406'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: converted

    script:
    """
    genome=\$(find staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
    samtools faidx "\$genome"
    for cram in staged/cram_*.cram; do
        sample=\$(basename "\${cram%.cram}")
        samtools view -@ ${task.cpus} -b -T "\$genome" "\$cram" > "staged/\${sample}.bam"
    done
    """

    stub:
    """
    mkdir -p staged
    for cram in staged/cram_*.cram; do
        sample=\$(basename "\${cram%.cram}")
        touch "staged/\${sample}.bam"
    done
    """
}
