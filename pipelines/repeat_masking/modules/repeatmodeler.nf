process REPEATMODELER_REPEATMODELER {
    tag "${meta.id}"
    label 'process_repeatmodeler'

    conda "bioconda::repeatmodeler=2.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.5--pl5321hdfd78af_1' :
        'biocontainers/repeatmodeler:2.0.5--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.families.fa"),        optional: true, emit: fasta
    tuple val(meta), path("*.families.stk"),        optional: true, emit: stk
    tuple val(meta), path("*-families.fa"),         optional: true, emit: fasta_alt
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    BuildDatabase -name ${prefix} -engine rmblast ${fasta}

    RepeatModeler \\
        -database ${prefix} \\
        -pa       ${task.cpus} \\
        -LTRStruct \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version 2>&1 | head -1 | sed 's/RepeatModeler version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.families.fa ${prefix}.families.stk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: 2.0.5
    END_VERSIONS
    """
}
