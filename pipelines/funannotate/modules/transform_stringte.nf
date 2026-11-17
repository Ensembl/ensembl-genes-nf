process TRANSFORM_STRINGTE {
    tag "$meta.assembly_accession:$meta.id"
    label 'process_medium'
    container params.gffread_container ?: '/nfs/production/flicek/ensembl/shared_data/singularity_images/cache/gffread_latest.sif'

    input:
    tuple val(meta), path(staged)

    output:
    tuple val(meta), path('staged'), emit: transformed
    path 'versions.yml', emit: versions

    script:
    def stringtie = meta.files.stringtie ? 'staged/stringtie.gtf' : ''
    """
    if [[ -n '${stringtie}' ]]; then
        genome=\$(find staged -maxdepth 1 -type f -name 'genome.fa*' | head -n1)
        [[ -s "\$genome" ]]
        [[ -s '${stringtie}' ]]
        gffread '${stringtie}' -g "\$genome" -w staged/stringtie_transcripts.fa
        [[ -s staged/stringtie_transcripts.fa ]]
    fi
    cat > versions.yml <<VERSIONS
"\${task.process}":
    gffread: \$(gffread --version 2>&1 | head -n1)
VERSIONS
    """

    stub:
    """
    ${meta.files.stringtie ? 'touch staged/stringtie_transcripts.fa' : ':'}
    touch versions.yml
    """
}
