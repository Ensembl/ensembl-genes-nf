process BIGBED_TO_BED {
    tag { "${meta.gca_accession}:${meta.aggregate_id}:${meta.source_entity_id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/ucsc-bigbedtobed:482--h0b57e2e_0'

    input:
    tuple val(meta), path(bigbed)

    output:
    tuple val(meta), path("${meta.safe_source_id}.junctions.bed"), emit: bed
    path "${meta.safe_source_id}.bigbed_to_bed.versions.yml", emit: versions

    script:
    """
    bigBedToBed ${bigbed} ${meta.safe_source_id}.junctions.bed
    test -s ${meta.safe_source_id}.junctions.bed || { echo 'Source splice BigBed has no rows' >&2; exit 1; }
    printf '"%s":\\n    bigBedToBed: 482\\n' '${task.process}' > ${meta.safe_source_id}.bigbed_to_bed.versions.yml
    """

    stub:
    """
    printf 'chr1\\t0\\t100\\tchr1:1-100:+\\t10\\t+\\t10\\t0\\tGT-AG\\n' > ${meta.safe_source_id}.junctions.bed
    printf '"%s":\\n    bigBedToBed: stub\\n' '${task.process}' > ${meta.safe_source_id}.bigbed_to_bed.versions.yml
    """
}
