process AGGREGATE_JUNCTIONS {
    tag { "${meta.gca_accession}:${meta.aggregate_id}:splice_junction" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(beds)
    path aggregator

    output:
    tuple val(meta), path("${meta.safe_id}.junctions.bed"), emit: bed
    path "${meta.safe_id}.splice.provenance.tsv", emit: provenance
    path "${meta.safe_id}.splice.versions.yml", emit: versions

    script:
    """
    python3 ${aggregator} --inputs ${beds.join(' ')} --output ${meta.safe_id}.junctions.bed
    printf 'aggregate_id\\tsource_count\\taggregation_method\\n${meta.aggregate_id}\\t${meta.source_count}\\t${meta.aggregation_method}\\n' > ${meta.safe_id}.splice.provenance.tsv
    printf '"%s":\\n    aggregate_splice_junctions: pipeline\\n' '${task.process}' > ${meta.safe_id}.splice.versions.yml
    """

    stub:
    """
    printf 'chr1\\t0\\t100\\tchr1:1-100:+\\t10\\t+\\t10\\t0\\tGT-AG\\n' > ${meta.safe_id}.junctions.bed
    printf 'aggregate_id\\tsource_count\\taggregation_method\\n${meta.aggregate_id}\\t${meta.source_count}\\t${meta.aggregation_method}\\n' > ${meta.safe_id}.splice.provenance.tsv
    printf '"%s":\\n    aggregate_splice_junctions: stub\\n' '${task.process}' > ${meta.safe_id}.splice.versions.yml
    """
}
