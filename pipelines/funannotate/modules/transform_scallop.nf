process TRANSFORM_SCALLOP {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container 'docker://getwilds/gffread:0.12.7'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: transformed
    path 'versions.yml', emit: versions

    script:
    def scallop = meta.files.scallop ? 'staged/scallop.gtf' : ''
    """
    if [[ -n '${scallop}' ]]; then
        genome=\$(find -L staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
        [[ -s "\$genome" ]]
        [[ -s '${scallop}' ]]
        gffread '${scallop}' -g "\$genome" -w staged/scallop_transcripts.fa
        [[ -s staged/scallop_transcripts.fa ]]
    fi
    cat > versions.yml <<VERSIONS
"${task.process}":
    gffread: \$(gffread --version 2>&1 | head -n1)
VERSIONS
    """

    stub:
    """
    ${meta.files.scallop ? 'touch staged/scallop_transcripts.fa' : ':'}
    touch versions.yml
    """
}
