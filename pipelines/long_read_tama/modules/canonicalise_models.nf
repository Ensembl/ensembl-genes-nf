process CANONICALISE_COMBINED_MODELS {
    tag "${backend}:${meta.id}:combined-models"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(backend), val(meta), path(tama_bed)

    output:
    tuple val(backend), val(meta), path("${backend}_combined_models.bed"), emit: bed
    path "${backend}_combined_models.sha256", emit: checksum
    path "${backend}_bed12_normalisation.tsv", emit: normalisation
    path 'versions.yml', emit: versions

    script:
    """
    normalise_bed12.py ${tama_bed} ${backend}_combined_models.bed ${backend}_bed12_normalisation.tsv
    sha256sum ${backend}_combined_models.bed > ${backend}_combined_models.sha256
    printf '"%s":\n    canonical_model_product: pipeline\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\t0\t4\t${backend}.stub.1\t0\t+\t0\t4\t0\t1\t4,\t0,\n' > ${backend}_combined_models.bed
    printf 'metric\tvalue\nmodels\t1\nzero_length_blocks_removed\t0\n' > ${backend}_bed12_normalisation.tsv
    sha256sum ${backend}_combined_models.bed > ${backend}_combined_models.sha256
    printf '"%s":\n    canonical_model_product: stub\n' '${task.process}' > versions.yml
    """
}
