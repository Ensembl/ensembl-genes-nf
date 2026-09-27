process PREPARE_MSA_INPUTS {
    label 'process_light'
    container 'oras://community.wave.seqera.io/library/python:3.12-slim'
    tag "${meta.id}"
    input:
    tuple val(meta), path(instances)
    output:
    tuple val(meta), path('msa_inputs.jsonl'), emit: manifest
    tuple val(meta), path('beds/*.bed'), emit: exon_beds, optional: true
    path 'versions.yml', emit: versions
    script:
    """
    prepare_msa_beds.py --instances ${instances} --outdir beds --manifest msa_inputs.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_msa_beds: 0.1.0
    END_VERSIONS
    """
    stub:
    """
    mkdir -p beds
    touch msa_inputs.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_msa_beds: 0.1.0
    END_VERSIONS
    """
}
