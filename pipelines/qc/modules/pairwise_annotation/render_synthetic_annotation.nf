process RENDER_SYNTHETIC_ANNOTATION {
    tag { "${meta.pair_id}_${source_label}" }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/synthetic_annotation" }, mode: 'copy', pattern: '*.synthetic.gff3'

    input:
    tuple val(meta), val(source_label), path(canonical_models)

    output:
    tuple val(meta), val(source_label), path("${source_label}.synthetic.gff3"), emit: annotation
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py render-synthetic --input '${canonical_models}' --source '${source_label}' --output '${source_label}.synthetic.gff3'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > '${source_label}.synthetic.gff3'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """
}
