process BLAST_BLASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::blast=2.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1' :
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta),  path(fasta)
    tuple val(meta2), path(db, stageAs: 'blast_db/*')
    val   out_ext

    output:
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    tuple val(meta), path("*.xml"), optional: true, emit: xml
    tuple val(meta), path("*.csv"), optional: true, emit: csv
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: meta.id
    def args    = task.ext.args   ?: ''
    def db_name = db[0].baseName
    """
    blastp \\
        -query    ${fasta} \\
        -db       blast_db/${db_name} \\
        -out      ${prefix}.blastp.${out_ext} \\
        -num_threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version 2>&1 | head -1 | sed 's/blastp: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.blastp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.15.0
    END_VERSIONS
    """
}
