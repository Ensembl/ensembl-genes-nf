process ENA_EXTRACT_BAM_HEADER {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path(file)

    output:
    tuple val(reference_key), val(meta), val(row), val(file_meta), path('source.bam'), path('bam.header'), emit: extracted
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    cp ${file} source.bam
    samtools view -H source.bam > bam.header
    printf 'ENA_EXTRACT_BAM_HEADER:\n  samtools: "%s"\n' "\$(samtools --version | awk 'NR==1 {print \$2}')" > versions.yml
    """

    stub:
    """
    touch source.bam bam.header
    printf 'ENA_EXTRACT_BAM_HEADER:\n  samtools: "stub"\n' > versions.yml
    """
}
