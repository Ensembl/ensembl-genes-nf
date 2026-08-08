process FIND_ADAPTERS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::biopython conda-forge::pandas"
    container "ghcr.io/lapti-ucc/riboseqorg-nf-python-pandas-sqlite:latest"

    input:
    tuple val(meta), path(raw_fastq)
    tuple val(meta2), path(fastqc_data)
    path adapter_list

    output:
    tuple val(meta), path("*_adapter_report.fa"), emit: adapter_report
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_adapters.py \\
        -i $fastqc_data \\
        -a ${adapter_list} \\
        -o "${prefix}_adapter_report.fa" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_adapters.py: \$(get_adapters.py --version 2>&1 | sed 's/get_adapters.py v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_adapter_report.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        get_adapters.py: 1.0.0
    END_VERSIONS
    """
}
