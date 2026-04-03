process REPEATMASKER_REPEATMASKER {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::repeatmasker=4.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_1' :
        'biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)
    path  lib        // repeat library file; pass [] to use -species instead

    output:
    tuple val(meta), path("*.masked"),      optional: true, emit: masked
    tuple val(meta), path("*.out"),                         emit: out
    tuple val(meta), path("*.gff"),         optional: true, emit: gff
    tuple val(meta), path("*.tbl"),         optional: true, emit: tbl
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def args     = task.ext.args   ?: ''
    def lib_flag = lib             ? "-lib ${lib}" : "-species '${params.species}'"
    """
    RepeatMasker \\
        ${lib_flag} \\
        -xsmall \\
        -gff \\
        -pa ${task.cpus} \\
        ${args} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v 2>&1 | head -1 | sed 's/RepeatMasker version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.fa.masked ${prefix}.fa.out ${prefix}.fa.gff ${prefix}.fa.tbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: 4.1.5
    END_VERSIONS
    """
}
