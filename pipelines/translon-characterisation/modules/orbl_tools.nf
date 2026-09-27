process ORBL_TOOLS {
    label 'process_medium'
    container params.orbl_container
    tag "${meta.id}"
    publishDir "${params.outdir}/06_constraint/orbl", mode: 'copy', pattern: '*.orbl.axis.jsonl'
    input:
    tuple val(meta), path(instances)
    output:
    tuple val(meta), path('*.orbl.axis.jsonl'), emit: axis
    path 'versions.yml', emit: versions
    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    prepare_orbl_input.py --instances ${instances} --output ${prefix}.orbl.input.tsv
    orbl.py ${params.orbl_alignment_set} --orblq --components ${prefix}.orbl.input.tsv > ${prefix}.orbl.raw.tsv
    parse_orbl_output.py --input ${prefix}.orbl.raw.tsv --output ${prefix}.orbl.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orbl_tools: \$(orbl.py --version 2>&1 | sed 's/.*Version: //')
        alignment_set: ${params.orbl_alignment_set}
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.orbl.axis.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orbl_tools: v0.9-beta
    END_VERSIONS
    """
}
