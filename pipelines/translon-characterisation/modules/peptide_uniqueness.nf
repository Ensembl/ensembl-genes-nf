process PEPTIDE_UNIQUENESS {
    label 'process_medium'
    container 'quay.io/biocontainers/mmseqs2:15-6f452--pl5321h6a68c12_0'
    tag "${meta.id}"
    publishDir "${params.outdir}/05_peptide_uniqueness", mode: 'copy', pattern: '*.uniqueness.jsonl'
    input:
    tuple val(meta), path(peptides)
    tuple val(proteome_meta), path(proteome)
    output:
    tuple val(meta), path('*.uniqueness.jsonl'), emit: results
    path 'versions.yml', emit: versions
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    peptide_uniqueness.py --peptides ${peptides} --proteome ${proteome} --output ${prefix}.uniqueness.jsonl ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.uniqueness.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: 15-6f452
    END_VERSIONS
    """
}
