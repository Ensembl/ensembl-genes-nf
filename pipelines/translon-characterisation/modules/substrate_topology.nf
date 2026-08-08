process SUBSTRATE_TOPOLOGY {
    label 'process_medium'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    publishDir "${params.outdir}/01_substrate", mode: 'copy', pattern: '*.instances.jsonl'

    input:
    tuple val(meta), path(intervals)
    tuple val(annotation_meta), path(gencode_gff3)
    tuple val(genome_meta), path(genome_fasta)

    output:
    tuple val(meta), path('*.instances.jsonl'), emit: instances
    tuple val(meta), path('*.unhosted.jsonl'), emit: unhosted
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    translon_characterise.py substrate --intervals ${intervals} --gff ${gencode_gff3} --genome ${genome_fasta} \\
        --instances ${prefix}.instances.jsonl --unhosted ${prefix}.unhosted.jsonl ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.instances.jsonl ${prefix}.unhosted.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
