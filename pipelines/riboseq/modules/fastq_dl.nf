process FASTQ_DL {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info/fastq_sources", mode: 'copy', pattern: '*.fastq_source.tsv'

    conda "bioconda::fastq-dl=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-dl:2.0.1--pyhdfd78af_0' :
        'biocontainers/fastq-dl:2.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(run), path(needs_processing)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "*.fastq_source.tsv", emit: source, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (params.exact_ena_fastq) {
        if (!(run ==~ /[SED]RR[0-9]+/)) {
            throw new IllegalArgumentException("Exact ENA acquisition requires a run accession: ${run}")
        }
        """
        fetch_ena_run.py ${run}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version | cut -d ' ' -f 2)
        END_VERSIONS
        """
    } else {
        def args = task.ext.args ?: ''
        """
        fastq-dl \\
            -a $run \\
            --cpus $task.cpus \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastq-dl: \$(fastq-dl --version | sed 's/fastq-dl //')
        END_VERSIONS
        """
    }

    stub:
    def files = params.exact_ena_fastq
        ? "${meta.id}.fastq.gz ${meta.id}.fastq_source.tsv"
        : "${meta.id}_1.fastq.gz ${meta.id}_2.fastq.gz"
    """
    touch ${files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: 2.0.1
    END_VERSIONS
    """
}
