process AUDIT_ANNOTATION_CONVERSION {
    tag { meta.pair_id }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/inputs" }, mode: 'copy', pattern: '*.tsv'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/synthetic_annotation" }, mode: 'copy', pattern: '*.id_map.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models)

    output:
    tuple val(meta), path('source_a_conversion_audit.tsv'), path('source_b_conversion_audit.tsv'), emit: audit
    tuple val(meta), path('source_a.id_map.tsv'), path('source_b.id_map.tsv'), emit: id_map
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py audit --models '${source_a_models}' --pair-id '${meta.pair_id}' --source source_a --source-format '${meta.format_a}' --output source_a_conversion_audit.tsv --id-map source_a.id_map.tsv
    pairwise_qc.py audit --models '${source_b_models}' --pair-id '${meta.pair_id}' --source source_b --source-format '${meta.format_b}' --output source_b_conversion_audit.tsv --id-map source_b.id_map.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    for source in source_a source_b; do
        printf 'pair_id\\tsource\\tsource_format\\tn_models\\tn_malformed\\tn_unique_genes\\tn_unique_loci\\tconversion_loss\\tstatus\\n' > "\${source}_conversion_audit.tsv"
        printf '%s\\t%s\\t%s\\t0\\t0\\t0\\t0\\t0\\tok\\n' '${meta.pair_id}' "\${source}" '${meta.format_a}' >> "\${source}_conversion_audit.tsv"
        printf 'pair_id\\tsource\\toriginal_gene_id\\toriginal_transcript_id\\tcanonical_transcript_id\\tgene_assignment_type\\tlocus_id\\n' > "\${source}.id_map.tsv"
    done
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
