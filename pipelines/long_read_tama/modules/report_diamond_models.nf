process REPORT_COMBINED_DIAMOND_MODELS {
    tag "${meta.id}:diamond-report"
    label 'process_medium'
    // This process runs the Python report generator; Diamond itself runs in
    // the upstream PREP_DIAMOND_DB/DIAMOND_BLASTP processes.
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(orf_manifest), path(diamond_hits)
    path diamond_db

    output:
    tuple val(meta), path('combined_diamond_hits.tsv'), emit: report
    path 'combined_diamond_summary.tsv', emit: summary
    path 'diamond_database_provenance.tsv', emit: provenance
    path 'versions.yml', emit: versions

    script:
    """
    report_diamond_models.py ${orf_manifest} ${diamond_hits} combined_diamond_hits.tsv combined_diamond_summary.tsv \\
        --classification-type ${params.diamond_classification_type}
    printf 'database_sha256\t%s\n' "\$(sha256sum ${diamond_db} | awk '{print \$1}')" > diamond_database_provenance.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: 2.1.24
    END_VERSIONS
    """

    stub:
    """
    printf 'model_id\tquery_peptide_id\tbest_subject_id\tpident\taligned_length\tquery_coverage\tevalue\tbitscore\tlegacy_classification\tlegacy_biotype_suffix\tstatus\n' > combined_diamond_hits.tsv
    printf 'model_id\tstatus\n' > combined_diamond_summary.tsv
    printf 'database_sha256\tstub\n' > diamond_database_provenance.tsv
    printf '"%s":\n    diamond: 2.1.24\n' '${task.process}' > versions.yml
    """
}
