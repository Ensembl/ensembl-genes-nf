process FASTQC {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::fastqc=0.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    path adapter_list

    output:
    tuple val(meta), path("*_fastqc.html", arity: '1'), emit: html
    tuple val(meta), path("*_fastqc.zip", arity: '1'), emit: zip
    tuple val(meta), path("*/*_fastqc_data.txt", arity: '1'), emit: txt
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fastqc
    fastqc \\
        $args \\
        --threads $task.cpus \\
        --outdir . \\
        --extract \\
        -q $fastq \\
        --adapters $adapter_list

    mv ${fastq.simpleName}_fastqc/fastqc_data.txt ${prefix}_fastqc/${prefix}_fastqc_data.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fastqc
    touch ${prefix}_fastqc.html
    touch ${prefix}_fastqc.zip
    touch ${prefix}_fastqc/${prefix}_fastqc_data.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}
