process AUDIT_ISOQUANT_OUTPUTS {
    tag "${meta.id}:${scope}:isoquant-audit"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(scope), path(products)

    output:
    tuple val(meta), val(scope), path('isoquant_manifest.tsv'), emit: manifest
    path 'isoquant_status.tsv', emit: status
    path 'versions.yml', emit: versions

    script:
    """
    audit_isoquant_outputs.py ${products} isoquant_manifest.tsv isoquant_status.tsv
    printf '"%s":\n    isoquant_output_audit: python\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'backend\tisoquant\nscope\t${scope}\ncomplete\ttrue\nstatus\tCOMPLETE\ntranscript_model_count\t1\nread_info_count\t1\nread_to_model_record_count\t1\n' > isoquant_manifest.tsv
    printf 'status\tcomplete\tmodels\tread_info\tread_to_model\nCOMPLETE\ttrue\t1\t1\t1\n' > isoquant_status.tsv
    printf '"%s":\n    isoquant_output_audit: stub\n' '${task.process}' > versions.yml
    """
}
