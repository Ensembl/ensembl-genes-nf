process AUDIT_LOCUS_ASSIGNMENT {
    tag { meta.pair_id }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/models" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models), path(assignment_summary, stageAs: 'assignment_input.tsv')

    output:
    tuple val(meta), path('locus_assignment_audit.tsv'), path('locus_assignment_summary.tsv'), emit: audit
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py audit-loci \\
      --source-a '${source_a_models}' --source-b '${source_b_models}' \\
      --pair-id '${meta.pair_id}' --assembly '${meta.assembly}' \\
      --reference '${meta.reference_annotation ?: ''}' \\
      --audit locus_assignment_audit.tsv --summary audit_summary.tsv
    cp assignment_input.tsv locus_assignment_summary.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tsource\\tassembly\\tlocus_id\\tgene_assignment_type\\tlocus_assignment_method\\tn_models\\ttranscript_ids\\n' > locus_assignment_audit.tsv
    cp assignment_input.tsv locus_assignment_summary.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """
}
