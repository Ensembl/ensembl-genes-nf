process VALIDATE_PAIR_MANIFEST {
    tag { meta.pair_id }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'

    input:
    tuple val(meta), val(source_a_ref), val(source_b_ref), val(assembly_report_ref), val(evidence_a_ref), val(evidence_b_ref)

    output:
    tuple val(meta), val(source_a_ref), val(source_b_ref), val(assembly_report_ref), val(evidence_a_ref), val(evidence_b_ref), path('validated_pair.tsv'), emit: validated
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py validate \\
      --pair-id '${meta.pair_id}' --assembly '${meta.assembly}' \\
      --source-a '${meta.source_a}' --source-b '${meta.source_b}' \\
      --source-a-path '${source_a_ref}' --source-b-path '${source_b_ref}' \\
      --source-a-format '${meta.format_a}' --source-b-format '${meta.format_b}' \\
      --source-a-assembly '${meta.source_a_assembly ?: ''}' --source-b-assembly '${meta.source_b_assembly ?: ''}' \\
      --output validated_pair.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tassembly\\tsource_a\\tsource_b\\tformat_a\\tformat_b\\n' > validated_pair.tsv
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' '${meta.pair_id}' '${meta.assembly}' '${meta.source_a}' '${meta.source_b}' '${meta.format_a}' '${meta.format_b}' >> validated_pair.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
