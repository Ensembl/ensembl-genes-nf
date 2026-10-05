process ISOQUANT_TO_BED12 {
    tag "${meta.id}:${scope}:isoquant-bed12"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(scope), path(gtf)

    output:
    tuple val(meta), val(scope), path('isoquant_models.bed'), emit: bed
    path 'versions.yml', emit: versions

    script:
    """
    test -s ${gtf} || { echo 'IsoQuant BED12 conversion received an empty GTF' >&2; exit 1; }
    gtf_to_bed12.py ${gtf} isoquant_models.bed '${meta.id}_${scope}'
    test -s isoquant_models.bed || { echo 'IsoQuant BED12 conversion produced no models' >&2; exit 1; }
    printf '"%s":\n    isoquant_to_bed12: python\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\t0\t4\tisoquant.stub.1\t0\t+\t0\t4\t0\t1\t4,\t0,\n' > isoquant_models.bed
    printf '"%s":\n    isoquant_to_bed12: stub\n' '${task.process}' > versions.yml
    """
}
