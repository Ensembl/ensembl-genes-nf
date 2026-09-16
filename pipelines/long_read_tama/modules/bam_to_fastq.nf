process BAM_TO_FASTQ {
    tag "${meta.id}:bam-to-fastq"
    label 'process_high'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(input_artifacts)

    output:
    tuple val(meta), path('canonical.fastq.gz'), emit: reads
    path 'versions.yml', emit: versions

    script:
    """
    input_bam=\$(find . -maxdepth 1 -type f -name '*.bam' -print -quit)
    test -n "\${input_bam}" || { echo "Expected CCS BAM for ${meta.id}" >&2; exit 1; }
    samtools fastq -F 0x900 "\${input_bam}" | gzip -c > canonical.fastq.gz
    test -s canonical.fastq.gz || { echo "BAM-to-FASTQ produced no reads for ${meta.id}" >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    printf '@stub/ccs\nACGT\n+\n!!!!\n' | gzip -c > canonical.fastq.gz
    printf '"%s":\n    samtools: 1.20\n' '${task.process}' > versions.yml
    """
}
