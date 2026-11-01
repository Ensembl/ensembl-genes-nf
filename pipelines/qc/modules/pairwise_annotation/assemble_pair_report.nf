process ASSEMBLE_PAIR_REPORT {
    tag { meta.pair_id }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/summary" }, mode: 'copy', pattern: 'pair_report.md'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/summary" }, mode: 'copy', pattern: 'versions.yml'

    input:
    tuple val(meta), path(structural_summary)

    output:
    tuple val(meta), path('pair_report.md'), emit: report
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py report --pair-id '${meta.pair_id}' --summary '${structural_summary}' --output pair_report.md
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf '# Pairwise annotation QC report\\n' > pair_report.md
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
