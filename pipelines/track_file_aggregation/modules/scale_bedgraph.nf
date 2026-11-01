process SCALE_BEDGRAPH {
    tag { "${meta.gca_accession}:${meta.aggregate_id}:${meta.coverage_strand}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(bedgraph)

    output:
    tuple val(meta), path("${meta.safe_id}.${meta.coverage_strand}.scaled.bedGraph"), emit: scaled
    path "${meta.safe_id}.${meta.coverage_strand}.scale.versions.yml", emit: versions

    script:
    """
    awk -v divisor=${meta.weight_divisor} 'BEGIN { OFS="\\t" } NF == 4 { \$4 = \$4 / divisor; print }' ${bedgraph} > ${meta.safe_id}.${meta.coverage_strand}.scaled.bedGraph
    test -s ${meta.safe_id}.${meta.coverage_strand}.scaled.bedGraph || { echo 'Scaled BigWig bedGraph is empty' >&2; exit 1; }
    printf '"%s":\\n    awk: repository-shell\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.scale.versions.yml
    """

    stub:
    """
    cp ${bedgraph} ${meta.safe_id}.${meta.coverage_strand}.scaled.bedGraph
    printf '"%s":\\n    awk: stub\\n' '${task.process}' > ${meta.safe_id}.${meta.coverage_strand}.scale.versions.yml
    """
}
