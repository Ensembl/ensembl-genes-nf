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
    test -s "${bed}" || { echo "Combined TAMA BED is missing or empty" >&2; exit 1; }
    awk 'NF != 12 {bad++} END {if (bad) {print "Invalid BED12 model rows: " bad > "/dev/stderr"; exit 1} print "models\\t" NR > "model_validation.tsv"}' "${bed}"
    printf '"%s":\\n    validator: awk\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'models\\t1\\n' > model_validation.tsv
    printf '"%s":\\n    validator: stub\\n' '${task.process}' > versions.yml
    """
}
