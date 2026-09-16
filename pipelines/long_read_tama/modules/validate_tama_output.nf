process VALIDATE_TAMA_OUTPUT {
    tag "${meta.id}:${shard}:tama-qc"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(shard), path(bed)

    output:
    tuple val(meta), val(shard), path('validated_tama.bed'), emit: bed
    path 'tama_validation.tsv', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    test -s "${bed}" || { echo "TAMA BED is missing or empty for ${meta.id}:${shard}" >&2; exit 1; }
    awk 'NF != 12 {bad++} END {if (bad) {print "Invalid TAMA BED12 rows: " bad > "/dev/stderr"; exit 1} print "run\\tshard\\tmodels\\n${meta.id}\\t${shard}\\t" NR > "tama_validation.tsv"}' "${bed}"
    cp -p "${bed}" validated_tama.bed
    printf '"%s":\n    tama_output_validation: pipeline\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\t0\t4\t${meta.id}.1\t0\t+\t0\t4\t0\t1\t4,\t0,\n' > validated_tama.bed
    printf 'run\tshard\tmodels\n${meta.id}\t${shard}\t1\n' > tama_validation.tsv
    printf '"%s":\n    tama_output_validation: stub\n' '${task.process}' > versions.yml
    """
}
