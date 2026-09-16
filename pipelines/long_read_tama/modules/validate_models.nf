process VALIDATE_LONG_READ_MODELS {
    tag 'models'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    path bed

    output:
    path 'model_validation.tsv', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    validate_tama_bed.py ${bed} model_validation.tsv
    printf '"%s":\\n    validator: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'models\\t1\\n' > model_validation.tsv
    printf '"%s":\\n    validator: stub\\n' '${task.process}' > versions.yml
    """
}
