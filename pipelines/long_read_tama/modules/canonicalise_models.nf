process CANONICALISE_COMBINED_MODELS {
    tag "${meta.id}:combined-models"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(tama_bed)

    output:
    tuple val(meta), path('combined_models.bed'), emit: bed
    path 'combined_models.sha256', emit: checksum
    path 'versions.yml', emit: versions

    script:
    """
    cp ${tama_bed} combined_models.bed
    sha256sum combined_models.bed > combined_models.sha256
    printf '"%s":\n    canonical_model_product: pipeline\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\t0\t4\tstub.1\t0\t+\t0\t4\t0\t1\t4,\t0,\n' > combined_models.bed
    sha256sum combined_models.bed > combined_models.sha256
    printf '"%s":\n    canonical_model_product: stub\n' '${task.process}' > versions.yml
    """
}
