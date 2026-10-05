process GTF_TO_BED12 {
    tag "${backend}:${meta.id}:${stage}:comparison-bed12"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(backend), val(stage), path(gtf)

    output:
    tuple val(meta), val(backend), val(stage), path('comparison_models.bed'), emit: bed
    path 'versions.yml', emit: versions

    script:
    """
    test -s ${gtf} || { echo 'GTF-to-BED12 conversion received an empty GTF' >&2; exit 1; }
    gtf_to_bed12.py ${gtf} comparison_models.bed '${meta.id}'
    test -s comparison_models.bed || { echo 'GTF-to-BED12 conversion produced no models' >&2; exit 1; }
    printf '"%s":\\n    gtf_to_bed12: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'chrStub\\t0\\t4\\t${meta.id}.stub.1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > comparison_models.bed
    printf '"%s":\\n    gtf_to_bed12: stub\\n' '${task.process}' > versions.yml
    """
}
