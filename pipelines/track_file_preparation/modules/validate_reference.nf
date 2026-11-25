process VALIDATE_REFERENCE {
    tag { "${assembly_release}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    path chrom_sizes
    val assembly_release
    path validator

    output:
    path 'reference_check.tsv', emit: checked
    path 'reference.versions.yml', emit: versions

    script:
    """
    python3 ${validator} --chrom-sizes ${chrom_sizes} --assembly-release '${assembly_release}'
    """

    stub:
    """
    printf 'assembly_release\\tchrom_sizes\\tstatus\\n${assembly_release}\\t${chrom_sizes}\\tcomplete\\n' > reference_check.tsv
    printf '"%s":\\n    python: stub\\n' '${task.process}' > reference.versions.yml
    """
}
