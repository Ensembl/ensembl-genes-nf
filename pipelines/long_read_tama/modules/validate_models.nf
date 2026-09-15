process VALIDATE_LONG_READ_MODELS {
    tag 'models'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11--he2b4eab_0'

    input:
    path bed

    output:
    path 'model_validation.tsv', emit: report

    script:
    """
    awk 'NF != 12 {bad++} END {if (bad) exit 1; print "models\\t" NR > "model_validation.tsv"}' ${bed}
    """

    stub:
    """
    printf 'models\\t1\\n' > model_validation.tsv
    """
}
