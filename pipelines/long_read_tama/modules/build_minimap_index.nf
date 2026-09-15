process BUILD_MINIMAP2_INDEX {
    tag "${reference.simpleName}"
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0'

    input:
    path reference

    output:
    path '*.mmi', emit: index
    path 'reference_checksums.tsv', emit: provenance
    path 'versions.yml', emit: versions

    script:
    def prefix = reference.baseName
    """
    minimap2 -d ${prefix}.mmi ${reference}
    sha256sum ${reference} ${prefix}.mmi > reference_checksums.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    def prefix = reference.baseName
    """
    touch ${prefix}.mmi
    sha256sum ${reference} ${prefix}.mmi > reference_checksums.tsv
    printf '"%s":\\n    minimap2: 2.28\\n' '${task.process}' > versions.yml
    """
}
