process ENA_CONVERT_TO_CRAM {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), val(row), val(file_meta), path(file)
    path reference_fasta
    path reference_fai

    output:
    tuple val(meta), val(row), val(file_meta), path('converted.cram'), emit: converted
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    samtools view -@ ${task.cpus} -C -T ${reference_fasta} -o converted.cram ${file}
    printf 'ENA_CONVERT_TO_CRAM:\n  samtools: "%s"\n' "\$(samtools --version | awk 'NR==1 {print \$2}')" > versions.yml
    """

    stub:
    """
    touch converted.cram
    printf 'ENA_CONVERT_TO_CRAM:\n  samtools: "stub"\n' > versions.yml
    """
}
