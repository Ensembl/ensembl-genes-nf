process AUDIT_NATIVE_SHARD_MODELS {
    tag "${backend}:${meta.id}:${shard}:model-count"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(backend), val(shard), path(native_model), val(native_format)

    output:
    tuple val(meta), val(backend), val(shard), path('audited_*'), val(native_format), emit: mergeable, optional: true
    tuple val(meta), val(backend), val(shard), path('shard_model_counts.tsv'), emit: counts
    path 'versions.yml', emit: versions

    script:
    def extension = native_format == 'bed12' ? 'bed' : 'gtf'
    def mergeableModel = "audited_${backend}_${meta.id}_${shard}.${extension}"
    """
    audit_native_shard_models.py '${native_model}' '${native_format}' '${backend}' '${meta.id}' '${shard}' '${mergeableModel}' shard_model_counts.tsv
    printf '"%s":\\n    native_shard_model_audit: python\\n' '${task.process}' > versions.yml
    """

    stub:
    def extension = native_format == 'bed12' ? 'bed' : 'gtf'
    def mergeableModel = "audited_${backend}_${meta.id}_${shard}.${extension}"
    """
    audit_native_shard_models.py '${native_model}' '${native_format}' '${backend}' '${meta.id}' '${shard}' '${mergeableModel}' shard_model_counts.tsv
    printf '"%s":\\n    native_shard_model_audit: python-stub\\n' '${task.process}' > versions.yml
    """
}

process SUMMARISE_NATIVE_SHARD_MODEL_COUNTS {
    tag 'native-contig-model-counts'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}/per_cohort/${params.cohort_id}/qc", mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }

    input:
    path(count_tables, stageAs: 'counts??/*')

    output:
    path 'native_contig_model_counts.tsv', emit: counts
    path 'versions.yml', emit: versions

    script:
    """
    summarise_native_shard_model_counts.py native_contig_model_counts.tsv
    printf '"%s":\\n    native_shard_model_count_summary: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'backend\\trun_accession\\tshard\\tnative_format\\tmodel_count\\tmerge_status\\n' > native_contig_model_counts.tsv
    printf '"%s":\\n    native_shard_model_count_summary: python-stub\\n' '${task.process}' > versions.yml
    """
}
