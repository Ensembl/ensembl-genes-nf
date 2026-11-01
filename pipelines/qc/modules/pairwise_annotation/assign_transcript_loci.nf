process ASSIGN_TRANSCRIPT_LOCI {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/models" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models), path(reference_annotation)

    output:
    tuple val(meta), path('source_a_canonical_models.tsv'), path('source_b_canonical_models.tsv'), path('source_a_locus_assignments.tsv'), path('source_b_locus_assignments.tsv'), path('locus_assignment_summary.tsv'), emit: loci
    path 'versions.yml', emit: versions

    script:
    def minimumOverlap = params.pairwise_locus_min_overlap_bp ?: 10
    def reciprocal = params.pairwise_locus_min_reciprocal_overlap ?: 0.10
    """
    set -euo pipefail
    pairwise_qc.py assign-loci \\
      --source-a '${source_a_models}' --source-b '${source_b_models}' \\
      --pair-id '${meta.pair_id}' --assembly '${meta.assembly}' \\
      --reference '${reference_annotation}' --reference-format '${meta.reference_format ?: 'auto'}' \\
      --minimum-overlap-bp '${minimumOverlap}' --reciprocal-overlap '${reciprocal}' \\
      --output-a source_a_canonical_models.tsv --output-b source_b_canonical_models.tsv \\
      --summary locus_assignment_summary.tsv
    cp source_a_canonical_models.tsv source_a_locus_assignments.tsv
    cp source_b_canonical_models.tsv source_b_locus_assignments.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """

    stub:
    """
    printf 'pair_id\\tsource\\tassembly\\tcontig\\tstrand\\ttranscript_id\\toriginal_transcript_id\\toriginal_gene_id\\ttranscript_start\\ttranscript_end\\texon_count\\texon_coordinates\\tintron_chain\\tfive_prime_end\\tthree_prime_end\\tsource_format\\tgene_assignment_type\\tlocus_id\\tlocus_assignment_method\\tmodel_status\\tutr_coordinates\\tcds_coordinates\\n' > source_a_canonical_models.tsv
    cp source_a_canonical_models.tsv source_b_canonical_models.tsv
    cp source_a_canonical_models.tsv source_a_locus_assignments.tsv
    cp source_b_canonical_models.tsv source_b_locus_assignments.tsv
    printf 'pair_id\\tassembly\\tsource\\tgene_assignment_type\\tlocus_assignment_method\\tn_models\\tn_loci\\ttotal_pair_loci\\treference_annotation\\treference_sha256\\tminimum_exonic_overlap_bp\\tminimum_reciprocal_exon_overlap\\n' > locus_assignment_summary.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.1.0
    END_VERSIONS
    """
}
