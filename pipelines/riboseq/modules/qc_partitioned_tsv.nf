/*
 * QC gate for partitioned unique-read TSV conversion.
 */

process QC_PARTITIONED_TSV {
    tag "${meta.id}"
    label 'process_low'
    cpus 1
    time '30.m'
    memory '1.GB'
    errorStrategy 'terminate'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"
    publishDir "${params.outdir}/matrix_partition_qc", mode: 'copy'

    input:
    tuple val(meta), path(stats_json)

    output:
    tuple val(meta), path("${meta.id}.partition_qc.json"), emit: qc_json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def qc_enabled = params.matrix_partition_qc_enabled == null ? true : params.matrix_partition_qc_enabled
    def min_total_records = task.ext.min_total_records ?: (qc_enabled ? (params.matrix_partition_qc_min_total_records ?: 1) : 0)
    def min_total_counts = task.ext.min_total_counts ?: (qc_enabled ? (params.matrix_partition_qc_min_total_counts ?: 1) : 0)
    def max_catch_all_unique_fraction = task.ext.max_catch_all_unique_fraction ?: (qc_enabled ? (params.matrix_partition_qc_max_n_unique_fraction ?: 0.05) : 1.0)
    def max_catch_all_count_fraction = task.ext.max_catch_all_count_fraction ?: (qc_enabled ? (params.matrix_partition_qc_max_n_count_fraction ?: 0.05) : 1.0)
    def action = task.ext.action ?: (qc_enabled ? (params.matrix_partition_qc_action ?: 'fail') : 'warn')
    """
    qc_partitioned_tsv.py \\
        ${stats_json} \\
        --output ${meta.id}.partition_qc.json \\
        --min-total-records ${min_total_records} \\
        --min-total-counts ${min_total_counts} \\
        --max-catch-all-unique-fraction ${max_catch_all_unique_fraction} \\
        --max-catch-all-count-fraction ${max_catch_all_count_fraction} \\
        --action ${action}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","passed":true,"action":"stub","failures":[],"warnings":[]}' > ${meta.id}.partition_qc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: unknown
    END_VERSIONS
    """
}
